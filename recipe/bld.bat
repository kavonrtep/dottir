@echo on
setlocal enableextensions

REM See build.sh — same rust-toolchain.toml workaround.
if exist rust-toolchain.toml move rust-toolchain.toml rust-toolchain.toml.bak

cargo build --release --workspace --locked
if errorlevel 1 exit /b 1

if not exist "%LIBRARY_BIN%" mkdir "%LIBRARY_BIN%"
copy /Y target\release\dottir.exe     "%LIBRARY_BIN%\dottir.exe"     || exit /b 1
copy /Y target\release\dottir-gui.exe "%LIBRARY_BIN%\dottir-gui.exe" || exit /b 1
