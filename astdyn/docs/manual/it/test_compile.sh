#!/bin/bash
for chap in 04 05 06 07 08 09 10 11 12 13 14 15; do
  echo "Testing up to chapter $chap..."
  # Create temp main with chapters up to $chap
  sed "/include{${chap}_/q" main_it.tex | sed 's/%\\include/\\include/g' > temp_main.tex
  echo "\\end{document}" >> temp_main.tex
  pdflatex -interaction=nonstopmode -output-directory=../../build/manual temp_main.tex > /dev/null 2>&1
  if [ $? -eq 0 ]; then
    echo "✓ Chapters 1-$chap OK"
  else
    echo "✗ Chapter $chap FAILED"
    break
  fi
done
rm -f temp_main.tex
