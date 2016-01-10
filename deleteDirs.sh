#!/bin/bash
echo "Enter the output directory you wish to DELETE: outputFiles/"
read outputDir
echo "Are you sure you want to DELETE the output directory $outputDir including ALL its contents?"
read yn
case $yn in
    [Yy]* ) rm -r outputFiles/$outputDir/;
	rm -r VMD/$outputDir/;;
    [Nn]* ) exit;;
    * ) echo "Please answer y/n.";;
esac
	