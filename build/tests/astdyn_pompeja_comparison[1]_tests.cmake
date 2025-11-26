add_test([=[PompejaComparisonTest.FullComparison]=]  [==[/Users/michelebigi/VisualStudio Code/GitHub/ITALOccultLibrary/build/tests/astdyn_pompeja_comparison]==] [==[--gtest_filter=PompejaComparisonTest.FullComparison]==] --gtest_also_run_disabled_tests)
set_tests_properties([=[PompejaComparisonTest.FullComparison]=]  PROPERTIES WORKING_DIRECTORY [==[/Users/michelebigi/VisualStudio Code/GitHub/ITALOccultLibrary/build/tests]==] SKIP_REGULAR_EXPRESSION [==[\[  SKIPPED \]]==])
set(  astdyn_pompeja_comparison_TESTS PompejaComparisonTest.FullComparison)
