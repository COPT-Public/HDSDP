if(NOT DEFINED SOLVER OR NOT DEFINED INSTANCE)
    message(FATAL_ERROR "SOLVER and INSTANCE are required")
endif()

execute_process(
    COMMAND "${SOLVER}" "${INSTANCE}"
    RESULT_VARIABLE result
    OUTPUT_VARIABLE output
    ERROR_VARIABLE error
    TIMEOUT 30
)

if(NOT "${result}" STREQUAL "0")
    message(FATAL_ERROR "sdpasolve failed (${result}):\n${output}\n${error}")
endif()

# The command-line wrapper can return zero even when it has not solved an SDP.
if(NOT output MATCHES "SDP Status:[ \t]*Primal dual optimal")
    message(FATAL_ERROR "sdpasolve did not report an optimal SDP solution:\n${output}\n${error}")
endif()

message(STATUS "sdpasolve solved ${INSTANCE} to primal-dual optimality")
