#include "BaseQualityHelper.h"
#include "StringArray.h"
#include "Parameters.h"
#include "IntArray.h"
#include "glfHandler.h"
#include "Error.h"

#include <math.h>
#include <time.h>

void StringToArray(const String & input, IntArray & values, int desired)
   {
   StringArray tokens;
   tokens.AddTokens(input, ',');

   values.Dimension(desired);
   values.Zero();

   if (tokens.Length())
      for (int i = 0; i < desired; i++)
         values[i] = tokens[i % tokens.Length()].AsInteger();
   }

double sq(double x) { return x * x; }

int main(int argc, char ** argv)
   {
   printf("glfMerge V1.0.2 -- Merge SNP calls based on .glf or .glz files\n");
   printf("(c) 2009 Goncalo Abecasis, Sebastian Zoellner, Yun Li\n\n");

   ParameterList pl;

   String   qualities = "30,30";
   String   minDepths = "1,1";
   String   maxDepths = "200,200";
   String   outfile = "merged.glf";

   bool verbose = false;

   IntArray qualityFilter;
   IntArray lowDepthFilter;
   IntArray highDepthFilter;

   BEGIN_LONG_PARAMETERS(longParameters)
      LONG_PARAMETER_GROUP("Map Quality Filter")
         LONG_STRINGPARAMETER("qualities", &qualities)
      LONG_PARAMETER_GROUP("Depth Filters")
         LONG_STRINGPARAMETER("minDepths", &minDepths)
         LONG_STRINGPARAMETER("maxDepths", &maxDepths)
      LONG_PARAMETER_GROUP("Output")
         LONG_STRINGPARAMETER("outfile", &outfile)
         LONG_PARAMETER("verbose", &verbose)
   END_LONG_PARAMETERS();

   pl.Add(new LongParameters("Options", longParameters));
   int argstart = pl.ReadWithTrailer(argc, argv) + 1;
   pl.Status();

   time_t t;
   time(&t);

   printf("Analysis started on %s\n", ctime(&t));
   fflush(stdout);

   int n = argc - argstart;
   argv += argstart;

   if (n == 0)
      error("No glf files listed at the end of command line\n");

   StringToArray(qualities, qualityFilter, n);
   StringToArray(minDepths, lowDepthFilter, n);
   StringToArray(maxDepths, highDepthFilter, n);

   glfHandler * glf = new glfHandler[n];

   for (int i = 0; i < n; i++)
      if (!glf[i].Open(argv[i]))
         error("Failed to open genotype likelihood file [%s]\n", argv[i]);

   printf("Calling genotypes for files ...\n");
   for (int i = 0; i < n; i++)
      if (glf[i].isOpen())
         printf("  %s\n", argv[i]);
   printf("\n");

   glfHandler output;
   output.Create(outfile);

   long long depth = 0, originalDepth = 0;
   long long sites = 0, originalSites = 0;

   while (glf[0].NextSection())
      {
      for (int i = 1; i < n; i++)
         {
         glf[i].NextSection();

         if (glf[0].maxPosition != glf[i].maxPosition || glf[0].label != glf[i].label)
            {
            error("Genotype files '%s' and '%s' are not compatible ...\n"
                "    File '%s' has section %s with %d entries ...\n"
                "    File '%s' section %s with %d entries ...\n",
                argv[0], argv[i],
                argv[0], (const char *) glf[0].label, glf[0].maxPosition,
                argv[i], (const char *) glf[i].label, glf[i].maxPosition);
            }
         }

      printf("Processing section %s with %d entries\n", (const char *) glf[0].label, glf[0].maxPosition);

      output.BeginSection(glf[0].label, glf[0].maxPosition);

      for (int i = 0; i < n; i++)
         glf[i].NextBaseEntry();

      int position = glf[0].position;
      char refBase = glf[0].data.refBase;
      for (int i = 1; i < n; i++)
         if (position > glf[i].position)
            {
            position = glf[i].position;
            refBase = glf[i].data.refBase;
            }

      while (position < glf[0].maxPosition)
         {
         output.data.recordType = 1;
         output.data.refBase = refBase;
         output.data.depth = 0;
         output.data.mapQuality = 0;
         output.data.minLLK = 0;

         for (int i = 0; i < 10; i++)
            output.data.lk[i] = 0;

         for (int i = 0; i < n; i++)
            if (glf[i].position == position &&
                glf[i].data.mapQuality >= qualityFilter[i] &&
                glf[i].data.depth >= (unsigned) lowDepthFilter[i] &&
               (glf[i].data.depth <= (unsigned) highDepthFilter[i] || highDepthFilter[i] == 0))
                {
                int deltaMinMap = output.data.lk[0] + glf[i].data.lk[0];

                for (int j = 1; j < 10; j++)
                  if (deltaMinMap > output.data.lk[j] + glf[i].data.lk[j])
                     deltaMinMap = output.data.lk[j] + glf[i].data.lk[j];

                output.data.minLLK += deltaMinMap + glf[i].data.minLLK;

                for (int j = 0; j < 10; j++)
                  if (output.data.lk[j] + glf[i].data.lk[j] - deltaMinMap < 255)
                     output.data.lk[j] += glf[i].data.lk[j] - deltaMinMap;
                  else
                     output.data.lk[j] = 255;

                output.data.mapQuality = (char) sqrt(
                  (  (sq(output.data.mapQuality) * output.data.depth) +
                     (sq(glf[i].data.mapQuality) * glf[i].data.depth) ) /
                     (output.data.depth + glf[i].data.depth + 1e-30) );
                output.data.depth += glf[i].data.depth;

               if (verbose)
                  {
                   printf("lk[%d] : { %d", i, glf[i].data.lk[0]);
                   for (int j = 1; j < 10; j++)
                     printf(", %d", glf[i].data.lk[j]);
                   printf("} [map: %d, depth: %d]\n", glf[i].data.mapQuality, glf[i].data.depth);
                   }
                }

         for (int i = 0; i < n; i++)
            if (glf[i].position == position)
               originalDepth += glf[i].data.depth;
         originalSites++;

         if (output.data.depth)
            {
            if (verbose)
               {
               printf("output : { %d", output.data.lk[0]);
               for (int j = 1; j < 10; j++)
                 printf(", %d", output.data.lk[j]);
               printf("} [map: %d, depth: %d]\n", output.data.mapQuality, output.data.depth);
               }

            output.WriteEntry(position);

            depth += output.data.depth;
            sites ++;
            }

         for (int i = 0; i < n; i++)
            if (glf[i].position == position)
               glf[i].NextBaseEntry();

         position = glf[0].position;
         refBase = glf[0].data.refBase;
         for (int i = 1; i < n; i++)
            if (position > glf[i].position)
               {
               position = glf[i].position;
               refBase = glf[i].data.refBase;
               }
         }

      output.EndSection();
      }

   printf("Combined file includes %ld bases covering %ld sites (%.1f coverage)\n",
           depth, sites, depth / (sites + 1e-30));
   printf("Original files include %ld bases covering %ld sites (%.1f coverage)\n",
           originalDepth, originalSites, originalDepth / (originalSites + 1e-30));

   for (int i = 0; i < n; i++)
      if (glf[i].isOpen())
         glf[i].Close();

   output.Close();
   }


