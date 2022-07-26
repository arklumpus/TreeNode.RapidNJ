@echo off

echo.
echo [104;97mDeleting previous build...[0m

for /f %%i in ('dir /a:d /b Release\RapidNJTests\*') do rd /s /q Release\RapidNJTests\%%i
del Release\RapidNJTests\* /s /f /q 1>nul

echo.
echo [104;97mCopying common resources...[0m

xcopy Resources\* Release\RapidNJTests\ /s /y /h >nul

echo.
echo Building with target [94mwin-x64[0m

cd RapidNJTestHost
dotnet publish -c Release /p:PublishProfile=Properties\PublishProfiles\win-x64.pubxml
cd ..

echo.
echo [104;97mCreating ZIP file...[0m

cd Release\RapidNJTests

move win-x64 RapidNJTests-win-x64
bash -c "zip -r RapidNJTests-win-x64.zip RapidNJTests-win-x64 >/dev/null"

for /f %%i in ('dir /a:d /b "RapidNJTests-win-x64"\*') do rd /s /q "RapidNJTests-win-x64"\%%i
del RapidNJTests-win-x64\* /s /f /q 1>nul
rmdir RapidNJTests-win-x64

cd ..\..

echo.
echo Building with target [94mlinux-x64[0m

cd RapidNJTestHost
dotnet publish -c Release /p:PublishProfile=Properties\PublishProfiles\linux-x64.pubxml
cd ..

echo.
echo [104;97mCreating tarball...[0m

cd Release\RapidNJTests

move linux-x64 RapidNJTests-linux-x64
bash -c "tar -czf RapidNJTests-linux-x64.tar.gz RapidNJTests-linux-x64"

for /f %%i in ('dir /a:d /b "RapidNJTests-linux-x64"\*') do rd /s /q "RapidNJTests-linux-x64"\%%i
del RapidNJTests-linux-x64\* /s /f /q 1>nul
rmdir RapidNJTests-linux-x64

cd ..\..

echo.
echo Building with target [94mmac-x64[0m

cd RapidNJTestHost
dotnet publish -c Release /p:PublishProfile=Properties\PublishProfiles\mac-x64.pubxml
cd ..

echo.
echo [104;97mCreating ZIP file...[0m

cd Release\RapidNJTests

move mac-x64 RapidNJTests-mac-x64
bash -c "zip -r RapidNJTests-mac-x64.zip RapidNJTests-mac-x64 >/dev/null"

for /f %%i in ('dir /a:d /b "RapidNJTests-mac-x64"\*') do rd /s /q "RapidNJTests-mac-x64"\%%i
del RapidNJTests-mac-x64\* /s /f /q 1>nul
rmdir RapidNJTests-mac-x64

cd ..\..

echo.
echo Building with target [94mmac-arm64[0m

cd RapidNJTestHost
dotnet publish -c Release /p:PublishProfile=Properties\PublishProfiles\mac-arm64.pubxml /p:PlatformTarget=arm64
cd ..

echo.
echo [104;97mCreating ZIP file...[0m

cd Release\RapidNJTests

move mac-arm64 RapidNJTests-mac-arm64
bash -c "zip -r RapidNJTests-mac-arm64.zip RapidNJTests-mac-arm64 >/dev/null"

for /f %%i in ('dir /a:d /b "RapidNJTests-mac-arm64"\*') do rd /s /q "RapidNJTests-mac-arm64"\%%i
del RapidNJTests-mac-arm64\* /s /f /q 1>nul
rmdir RapidNJTests-mac-arm64

cd ..\..

echo.
echo Building with target [94mlinux-arm64[0m

cd RapidNJTestHost
dotnet publish -c Release /p:PublishProfile=Properties\PublishProfiles\linux-arm64.pubxml /p:PlatformTarget=arm64
cd ..

echo.
echo [104;97mCreating tarball...[0m

cd Release\RapidNJTests

move linux-arm64 RapidNJTests-linux-arm64
bash -c "tar -czf RapidNJTests-linux-arm64.tar.gz RapidNJTests-linux-arm64"

for /f %%i in ('dir /a:d /b "RapidNJTests-linux-arm64"\*') do rd /s /q "RapidNJTests-linux-arm64"\%%i
del RapidNJTests-linux-arm64\* /s /f /q 1>nul
rmdir RapidNJTests-linux-arm64

cd ..\..

echo.
echo [94mAll done![0m