.PHONY: task clean


task:
	module load r/4.4.0; Rscript R/task.R

clean:
	rm -f results/*.rds results/*.png figs/*.html
