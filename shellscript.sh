#!/bin/sh
echo "Executing carbonmap.R script in the background"
Rscript --vanilla carbonmap.R > carbonmap.log 2>&1 &
echo "Check the progress with command 'tail -f carbonmap.log'"
echo "Check the processor usage with command 'top'"
## End of script