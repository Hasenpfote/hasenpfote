
#
if(WIN32)
	if(MSVC)
		#[[
		if(CMAKE_CONFIGURATION_TYPES)
			set(CMAKE_CONFIGURATION_TYPES "Debug;Release" CACHE STRING "Multi config types" FORCE)
		endif()
		]]
	endif()
elseif(UNIX)
	if(NOT CMAKE_BUILD_TYPE)
		set(CMAKE_BUILD_TYPE Release CACHE STRING "Force a Release build." FORCE)
	endif()
	set(CMAKE_CXX_FLAGS "-std=c++14 -Wall -Wextra -pedantic -Wcast-align -Wcast-qual -Wconversion -Wdisabled-optimization -Wendif-labels -Wfloat-equal -Winit-self -Winline -Wlogical-op -Wmissing-include-dirs -Wnon-virtual-dtor -Wold-style-cast -Woverloaded-virtual -Wpacked -Wpointer-arith -Wredundant-decls -Wshadow -Wsign-promo -Wswitch-default -Wswitch-enum -Wunsafe-loop-optimizations -Wvariadic-macros -Wwrite-strings ")
	set(CMAKE_CXX_FLAGS_DEBUG "-g3 -O0 -pg")
	set(CMAKE_CXX_FLAGS_RELEASE "-O2 -s -DNDEBUG -march=native")
	set(CMAKE_CXX_FLAGS_MINSIZEREL "-Os -s -DNDEBUG -march=native")
	set(CMAKE_CXX_FLAGS_RELWITHDEBINFO "-g3 -Og -pg")
endif()

#
add_definitions(-D_UNICODE -DUNICODE)
