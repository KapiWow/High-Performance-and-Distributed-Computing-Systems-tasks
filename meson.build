project('atata', 'cpp')

cpp = meson.get_compiler('cpp')

srcs = [
    'task1',
    'task2_1',
    'task2_2'
]

srcs_openmp = [
    'task2_3',
]

lib = [
]

cpp_arg = [
    '-O4',
    '-march=native',
    '-ffast-math'
]

openmp_arg = [
	'-lpthread',
    '-fopenmp',
]

linker_arg = [
    '-fopenmp',
]

foreach src : srcs
	src_name= src + '.cc'
	executable(
	    src, 
	    sources: src_name,
	    dependencies : lib,
	    cpp_args : cpp_arg,
	)
endforeach

foreach src : srcs_openmp
	src_name= src + '.cc'
	executable(
	    src, 
	    sources: src_name,
	    dependencies : lib,
	    cpp_args : [cpp_arg + openmp_arg],
        link_args : linker_arg,
	)
endforeach

