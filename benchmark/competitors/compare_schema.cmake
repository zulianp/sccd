# The competitor harness and sccd_bench must print the same CSV schema.
#
# A competitor row exists to sit beside an SCCD row in one file. If the two
# disagree on the columns, the report reads one of them with the wrong names and
# says nothing about it -- which is exactly how every toi_* column was once
# dropped unnamed. Each binary prints its own header; this fails the build's test
# suite when they drift apart.

execute_process(COMMAND "${SCCD_BENCH}" --header
                OUTPUT_VARIABLE sccd_header OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE sccd_rc)
execute_process(COMMAND "${COMPETITOR_BENCH}" --header
                OUTPUT_VARIABLE competitor_header OUTPUT_STRIP_TRAILING_WHITESPACE
                RESULT_VARIABLE competitor_rc)

if(NOT sccd_rc EQUAL 0)
  message(FATAL_ERROR "sccd_bench --header failed (${sccd_rc})")
endif()
if(NOT competitor_rc EQUAL 0)
  message(FATAL_ERROR "competitor --header failed (${competitor_rc})")
endif()

if(NOT sccd_header STREQUAL competitor_header)
  message(FATAL_ERROR
    "CSV schemas disagree.\n"
    "  sccd_bench: ${sccd_header}\n"
    "  competitor: ${competitor_header}\n"
    "Update kCsvHeader in benchmark/competitors/common/competitor_bench.hpp.")
endif()

message(STATUS "schemas agree: ${sccd_header}")
