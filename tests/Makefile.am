# End-to-end workflow tests: run each workflow on small inputs and compare
# every stage's output with reference outputs generated from master.
TESTS = ir.test raman.test sfg.test power.test
TEST_EXTENSIONS = .test
TEST_LOG_COMPILER = $(SHELL)
AM_TESTS_ENVIRONMENT = top_builddir='$(abs_top_builddir)'; export top_builddir;

EXTRA_DIST = \
    $(TESTS) \
    README.md \
    common.sh \
    compare.sh \
    generate-reference.sh \
    data \
    reference
