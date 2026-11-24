# Vendored Autoconf macros

`ax_cxx_compile_stdcxx.m4` is copied without modification from the
[GNU Autoconf Archive at revision e5369595](https://github.com/autoconf-archive/autoconf-archive/blob/e5369595076dbe5b8d0ab2d255735dc31f9060e9/m4/ax_cxx_compile_stdcxx.m4).
Its copyright and redistribution terms are retained in the file.

This revision includes the standard-mode selection fixes following the
2024.10.16 release and predates a test comment containing shell backticks,
which Autoconf expands as command substitution in its generated heredoc.
Keep the macro pinned and test generated configure scripts when updating it.

`ax_check_link_flag.m4` (serial 6) is copied without modification from the
GNU Autoconf Archive release 2024.10.16, as packaged by Ubuntu's
`autoconf-archive` 20241016-3. `configure` uses it to check
`--enable-sanitizers`.

`pkg.m4` (serial 13) is copied without modification from pkgconf 2.5.1, as
packaged by Ubuntu's `pkgconf` 2.5.1-4. Vendoring it lets `autogen.sh` run
without pkg-config installed; pkg-config itself is still optional when running
`configure`. Its default action when pkg-config is missing is to abort, so
`configure.ac` passes an explicit empty action.
