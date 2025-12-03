.PHONY: all clean build run

# Default target
all: report/03_report.pdf plot/Figure4.png

plot/Figure4.png: src/01_exploration.r data/poker-hand-testing.data data/poker-hand-training-true.data
	mkdir -p plot
	Rscript src/01_exploration.r

saved_model/latest_model.model: src/02_xgboost_model.r data/poker-hand-testing.data data/poker-hand-training-true.data
	mkdir -p saved_model
	Rscript src/02_xgboost_model.r

report/03_report.pdf: report/03_report.Rmd saved_model/latest_model.model
	Rscript -e "rmarkdown::render('report/03_report.Rmd', output_format = 'pdf_document')"

clean:
	rm -f report/03_report.pdf
	rm -f plot/*.png
	rm -f saved_model/latest_model.model

build:
	docker build . -t zehui

run:
	docker run --rm -v $(PWD):/home/rstudio/project -w /home/rstudio/project zehui make all
