project('atata', 'cpp')

cpp = meson.get_compiler('cpp')

srcs = [
    'task1',
]

lib = [
]

cpp_arg = [
    '-O4',
    '-march=native',
#	'-g',
#	'-masm=intel',
#	'-lpthread',
#	'-ffast-math',
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

