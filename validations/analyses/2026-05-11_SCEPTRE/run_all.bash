RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI1_1_sceptre_object_low_moi_merged.rds', \
              prefix='MOI1_1', integration_strategy='union', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI1_1_union.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI1_5_sceptre_object_low_moi_merged.rds', \
              prefix='MOI1_5', integration_strategy='union', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI1_5_union.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI11_3_sceptre_object_high_moi_merged.rds', \
              prefix='MOI11_3', integration_strategy='union', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI11_3_union.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI13_sceptre_object_high_moi_merged.rds', \
              prefix='MOI13', integration_strategy='union', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI13_union.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI1_1_sceptre_object_low_moi_merged.rds', \
              prefix='MOI1_1', integration_strategy='singleton', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI1_1_singleton.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI1_5_sceptre_object_low_moi_merged.rds', \
              prefix='MOI1_5', integration_strategy='singleton', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI1_5_singleton.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI11_3_sceptre_object_high_moi_merged.rds', \
              prefix='MOI11_3', integration_strategy='singleton', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI11_3_singleton.html')" 2>&1 | tail -30

RSTUDIO_PANDOC='/Applications/RStudio.app/Contents/Resources/app/quarto/bin/tools/aarch64' \
/usr/local/bin/Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/MOI13_sceptre_object_high_moi_merged.rds', \
              prefix='MOI13', integration_strategy='singleton', n_processors=4), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_MOI13_singleton.html')" 2>&1 | tail -30
  
Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/low_moi_merged_sceptre_object_low_moi_merged.rds', \
              prefix='low_moi_merged', integration_strategy='union', n_processors=2), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_low_moi_merged_union.html')" 2>&1 | tail -30

Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/low_moi_merged_sceptre_object_low_moi_merged.rds', \
              prefix='low_moi_merged', integration_strategy='singleton', n_processors=2), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_low_moi_merged_singleton.html')" 2>&1 | tail -30

Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/high_moi_merged_sceptre_object_high_moi_merged.rds', \
              prefix='high_moi_merged', integration_strategy='union', n_processors=2), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_high_moi_merged_union.html')" 2>&1 | tail -30

Rscript -e "rmarkdown::render('analyses/2026-05-11_SCEPTRE/02_run_sceptre_pipeline.Rmd', \
  params=list(sceptre_object_path='results/2026-05-11_SCEPTRE/objects/high_moi_merged_sceptre_object_high_moi_merged.rds', \
              prefix='high_moi_merged', integration_strategy='singleton', n_processors=2), \
  output_format='html_document', output_dir='results/2026-05-11_SCEPTRE', \
  output_file='02_high_moi_merged_singleton.html')" 2>&1 | tail -30