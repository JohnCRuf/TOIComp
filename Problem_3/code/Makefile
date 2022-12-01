SHELL=bash
include ../../shell_functions.make

outputs_b = ../output/problem_3_part_b_value_comparison_delta_0p3.png ../output/problem_3_part_b_table.png
outputs_c = ../output/problem_3_part_c_value_comparison_delta_0p3.png
outputs =  $(outputs_b) $(outputs_c)
all: $(outputs)
b: $(outputs_b)
c: $(outputs_c)

$(outputs_b): part_b.jl | ../output problem_3_fns.jl
	$(JULIA) $<
$(outputs_c): part_c.jl | ../output problem_3_fns.jl
	$(JULIA) $<

include ../../generic.make