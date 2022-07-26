cmake_minimum_required (VERSION 3.8)

add_library (rapidNJWrapper SHARED "rapidNJWrapper.cpp" "rapidNJWrapper.h" )

target_include_directories(rapidNJWrapper PRIVATE "include/all")
target_include_directories(rapidNJWrapper PRIVATE ".")

if (WIN32)
	target_include_directories(rapidNJWrapper PRIVATE "include/win-all")
endif()

set(TARGET "AMD64" CACHE STRING "Target build architecture for Windows.")

if (WIN32)
	if (${TARGET} MATCHES "AMD64")
		target_link_libraries(rapidNJWrapper "${CMAKE_CURRENT_SOURCE_DIR}/lib/win-x64/rapidNJ.lib")
		target_link_libraries(rapidNJWrapper "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/rapidNJ/lib/x64/pthreadVC2S.lib")
	endif()
	if (${TARGET} MATCHES "x86")
		target_link_libraries(rapidNJWrapper "${CMAKE_CURRENT_SOURCE_DIR}/lib/win-x64/rapidNJ.lib")
		target_link_libraries(rapidNJWrapper "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/rapidNJ/lib/x64/pthreadVC2S.lib")
	endif()
	if (${TARGET} MATCHES "ARM64")
		target_link_libraries(rapidNJWrapper "${CMAKE_CURRENT_SOURCE_DIR}/lib/win-x64/rapidNJ.lib")
		target_link_libraries(rapidNJWrapper "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/rapidNJ/lib/x64/pthreadVC2S.lib")
	endif()
endif()

if(UNIX AND NOT APPLE)
	if (${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
		target_link_libraries(rapidNJWrapper "${CMAKE_CURRENT_SOURCE_DIR}/lib/linux-x64/rapidNJ.a")
	endif()
	
	if (${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "aarch64")
		target_link_libraries(rapidNJWrapper "${CMAKE_CURRENT_SOURCE_DIR}/lib/win-x64/rapidNJ.lib")
		target_link_libraries(rapidNJWrapper "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/rapidNJ/lib/x64/pthreadVC2S.lib")
	endif()
endif()

if (APPLE)
	if (${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "x86_64")
		target_link_libraries(rapidNJWrapper "${CMAKE_CURRENT_SOURCE_DIR}/lib/mac-x64/rapidNJ.a")
	endif()
	
	if (${CMAKE_HOST_SYSTEM_PROCESSOR} MATCHES "arm64")
		target_link_libraries(rapidNJWrapper "${CMAKE_CURRENT_SOURCE_DIR}/lib/win-x64/rapidNJ.lib")
		target_link_libraries(rapidNJWrapper "${CMAKE_CURRENT_SOURCE_DIR}/thirdparty/rapidNJ/lib/x64/pthreadVC2S.lib")
	endif()
endif()