#include <simdjson.h>

#include <cstdint>
#include <cstdlib>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <limits>
#include <vector>

namespace fs = std::filesystem;

namespace {

bool outputs_are_current(const fs::path& input, const fs::path& c0_path, const fs::path& c1_path)
{
    std::error_code ec;
    if (!fs::exists(c0_path, ec) || !fs::exists(c1_path, ec)) {
        return false;
    }
    const auto input_time = fs::last_write_time(input, ec);
    if (ec) {
        return false;
    }
    const auto c0_time = fs::last_write_time(c0_path, ec);
    if (ec || c0_time < input_time) {
        return false;
    }
    const auto c1_time = fs::last_write_time(c1_path, ec);
    return !ec && c1_time >= input_time;
}

template <typename T>
bool write_raw(const fs::path& path, const std::vector<T>& values)
{
    std::ofstream out(path, std::ios::binary | std::ios::trunc);
    if (!out) {
        std::cerr << "error: failed to open " << path << " for writing\n";
        return false;
    }

    const auto bytes = static_cast<std::streamsize>(values.size() * sizeof(T));
    out.write(reinterpret_cast<const char*>(values.data()), bytes);
    if (!out) {
        std::cerr << "error: failed to write " << path << "\n";
        return false;
    }

    return true;
}

bool parse_pair(simdjson::dom::element row, std::int32_t& c0, std::int32_t& c1)
{
    simdjson::dom::array pair;
    auto error = row.get_array().get(pair);
    if (error || pair.size() != 2) {
        return false;
    }

    std::int64_t first = 0;
    std::int64_t second = 0;
    error = pair.at(0).get_int64().get(first);
    if (error) {
        return false;
    }
    error = pair.at(1).get_int64().get(second);
    if (error) {
        return false;
    }

    constexpr auto int32_min = static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::min());
    constexpr auto int32_max = static_cast<std::int64_t>(std::numeric_limits<std::int32_t>::max());
    if (first < int32_min || first > int32_max || second < int32_min || second > int32_max) {
        return false;
    }

    c0 = static_cast<std::int32_t>(first);
    c1 = static_cast<std::int32_t>(second);
    return true;
}

bool convert_file(simdjson::dom::parser& parser, const fs::path& input)
{
    const fs::path output_dir = input.parent_path() / input.stem();
    const fs::path c0_path = output_dir / "c0.int32";
    const fs::path c1_path = output_dir / "c1.int32";

    if (outputs_are_current(input, c0_path, c1_path)) {
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

    std::vector<std::int32_t> c0;
    std::vector<std::int32_t> c1;
    c0.reserve(rows.size());
    c1.reserve(rows.size());

    for (simdjson::dom::element row : rows) {
        std::int32_t first = 0;
        std::int32_t second = 0;
        if (!parse_pair(row, first, second)) {
            std::cerr << "error: expected an int32 pair in " << input << "\n";
            return false;
        }
        c0.push_back(first);
        c1.push_back(second);
    }

    std::error_code ec;
    fs::create_directories(output_dir, ec);
    if (ec) {
        std::cerr << "error: failed to create " << output_dir << ": " << ec.message() << "\n";
        return false;
    }

    return write_raw(c0_path, c0) && write_raw(c1_path, c1);
}

} // namespace

int main(int argc, char** argv)
{
    if (argc < 2) {
        std::cerr << "usage: " << argv[0] << " <boxes.json> [<boxes.json> ...]\n";
        return EXIT_FAILURE;
    }

    simdjson::dom::parser parser;
    bool ok = true;
    for (int i = 1; i < argc; ++i) {
        ok = convert_file(parser, argv[i]) && ok;
    }

    return ok ? EXIT_SUCCESS : EXIT_FAILURE;
}
