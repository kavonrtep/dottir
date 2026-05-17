@echo off
REM Move rust-toolchain.toml aside so cargo uses the conda-provided
REM rust rather than asking rustup for the pinned channel.
if exist rust-toolchain.toml move rust-toolchain.toml rust-toolchain.toml.bak

cargo build --release --workspace --locked
if errorlevel 1 exit /b 1

copy /Y target\release\dottir.exe     "%LIBRARY_BIN%\dottir.exe"
if errorlevel 1 exit /b 1
copy /Y target\release\dottir-gui.exe "%LIBRARY_BIN%\dottir-gui.exe"
if errorlevel 1 exit /b 1
