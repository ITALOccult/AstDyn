add_test([=[ResidualDebug.SimpleTest]=]  [==[/Users/michelebigi/VisualStudio Code/GitHub/ITALOccultLibrary/astdyn/tests/astdyn_residual_debug]==] [==[--gtest_filter=ResidualDebug.SimpleTest]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[ResidualDebug.SimpleTest]=]  PROPERTIES WORKING_DIRECTORY [==[/Users/michelebigi/VisualStudio Code/GitHub/ITALOccultLibrary/astdyn/tests]==] SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set(  astdyn_residual_debug_TESTS ResidualDebug.SimpleTest)
