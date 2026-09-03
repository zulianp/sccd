#include <simdjson.h>

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

namespace {

bool output_is_current(const fs::path& input, const fs::path& output)
{
    std::error_code ec;
    if (!fs::exists(output, ec)) {
        return false;
    }

    const auto input_time = fs::last_write_time(input, ec);
    if (ec) {
        return false;
    }
    const auto output_time = fs::last_write_time(output, ec);
    return !ec && output_time >= input_time;
}

bool write_raw(const fs::path& path, const std::vector<std::uint8_t>& values)
{
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        std::cerr << "error: failed to open " << path << " for writing\n";
        return false;
    }

    out.write(reinterpret_cast<const char*>(values.data()), static_cast<std::streamsize>(values.size()));
    if (!out) {
        std::cerr << "error: failed to write " << path << "\n";
        return false;
    }

    return true;
}

fs::path output_dir_for(const fs::path& input)
{
    const std::string suffix = "_mma_bool";
    std::string key = input.stem().string();
    if (key.size() >= suffix.size() && key.compare(key.size() - suffix.size(), suffix.size(), suffix) == 0) {
        key.resize(key.size() - suffix.size());
    }
    return input.parent_path() / key;
}

bool convert_file(simdjson::dom::parser& parser, const fs::path& input)
{
    const fs::path output_dir = output_dir_for(input);
    const fs::path output_path = output_dir / "mma_bool.uint8";
    if (output_is_current(input, output_path)) {
        return true;
    }

    simdjson::dom::element doc;
    auto error = parser.load(input.string()).get(doc);
    if (error) {
        std::cerr << "error: failed to parse " << input << ": " << error << "\n";
        return false;
    }

    simdjson::dom::array rows;
    error = doc.get_array().get(rows);
    if (error) {
        std::cerr << "error: expected top-level array in " << input << "\n";
        return false;
    }

    std::vector<std::uint8_t> values;
    values.reserve(rows.size());
    for (simdjson::dom::element row : rows) {
        bool value = false;
        error = row.get_bool().get(value);
        if (error) {
            std::cerr << "error: expected a boolean in " << input << "\n";
            return false;
        }
        values.push_back(value ? 1 : 0);
    }

    std::error_code ec;
    fs::create_directories(output_dir, ec);
    if (ec) {
        std::cerr << "error: failed to create " << output_dir << ": " << ec.message() << "\n";
        return false;
    }

    return write_raw(output_path, values);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <mma_bool.json> [<mma_bool.json> ...]\n";
        return EXIT_FAILURE;
    }

    simdjson::dom::parser parser;
    bool ok = true;
    for (int i = 1; i < argc; ++i) {
        ok = convert_file(parser, argv[i]) && ok;
    }

    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
