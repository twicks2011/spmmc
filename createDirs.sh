#!/bin/bash
echo "Enter output directory you wish to CREATE: outputFiles/"
read outputDir

mkdir outputFiles/$outputDir/
mkdir VMD/$outputDir/

mkdir outputFiles/$outputDir/occupancies/
mkdir outputFiles/$outputDir/restartInfo/
