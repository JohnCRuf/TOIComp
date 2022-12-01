SHELL=bash
include ../../shell_functions.make

outputs_a = ../output/part_a_value_comparison.png ../output/part_a_policy_comparison.png
outputs_b = ../output/part_b_value_comparison_delta_0p3.png ../output/part_b_table.png
outputs_c = ../output/part_c_computation_time.png ../output/part_c_error.png
outputs = $(outputs_a) $(outputs_b) $(outputs_c)
all: $(outputs)

$(outputs_a): part_a.jl | ../output problem_2_fns.jl
	$(JULIA) $<
$(outputs_b): part_b.jl | ../output problem_2_fns.jl
	$(JULIA) $<
$(outputs_c): part_c.jl | ../output problem_2_fns.jl
	$(JULIA) $<

include ../../generic.make