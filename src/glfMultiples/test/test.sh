ERROR=false
../../bin/glfMultiples --minMapQuality 0 --minDepth 1 --maxDepth 10000000 --uniformPrior --glfAliases testFiles/glfAlias.txt -b results/output.vcf testFiles/NA06984.20.1.1000000.glf testFiles/NA06985.20.1.1000000.glf testFiles/NA06986.20.1.1000000.glf testFiles/NA06989.20.1.1000000.glf testFiles/NA06994.20.1.1000000.glf testFiles/NA07000.20.1.1000000.glf testFiles/NA07037.20.1.1000000.glf > results/output.txt 2>&1 \
&& diff -I "##filedate=.*" results/output.vcf expected/expected.vcf \
&& diff -I "^Analysis .*ed on [0-9:]*" results/output.txt expected/expected.txt

if [ $? -ne 0 ]
then
    ERROR=true
fi

if($ERROR == true)
then
  exit 1
fi
