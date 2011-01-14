////////////////////////////////////////////////////////////////////// 
// pdf/PDFpage.h 
// (c) 2000-2010 Goncalo Abecasis
// 
// This file is distributed as part of the Goncalo source code package   
// and may not be redistributed in any form, without prior written    
// permission from the author. Permission is granted for you to       
// modify this file for your own personal use, but modified versions  
// must retain this copyright notice and must not be distributed.     
// 
// Permission is granted for you to use this file to compile Goncalo.    
// 
// All computer programs have bugs. Use this file at your own risk.   
// 
// Wednesday June 16, 2010
// 
 
#ifndef __PDFPAGE_H__
#define __PDFPAGE_H__

#include "IntArray.h"
#include "PDFfont.h"

class PDF;

enum PDFPageSize   {psLetter, psLetterR, psA4, psA4R};
enum PDFLineCap    {lcButt = 0, lcRound = 1, lcSquare = 2};
enum PDFLineStyle  {lsSolid = 0, lsDashed = 1, lsDotted = 2};
enum PDFTextAlignH {taLeft, taRight, taCenter};
enum PDFTextAlignV {taAbove, taBelow, taMiddle};

class PDFPage
   {
   public:
      int tree_index;

      PDFTextAlignH hTextAlignment;
      PDFTextAlignV vTextAlignment;

      PDFPage(PDF & parent);

      void OpenPage();
      void ClosePage();
      void WritePageTree();

      void SetSize(PDFPageSize size);
      void SetSize(int width_in_points, int height_in_points);

      double GetHeight();
      double GetWidth();

      void SetLineColor(double red, double green, double blue);
      void SetFillColor(double red, double green, double blue);
      void SetLineCMYK(double cyan, double magenta, double yellow, double black);
      void SetFillCMYK(double cyan, double magenta, double yellow, double black);
      void SetLineGray(double gray);
      void SetFillGray(double gray);

      void DrawLine(double x0, double y0, double x1, double y1);

      void DrawRectangle(double x0, double y0, double x1, double y1);
      void FillRectangle(double x0, double y0, double x1, double y1);
      void Rectangle(double x0, double y0, double x1, double y1);

      void DrawPolygon(double * x, double * y, int points);
      void FillPolygon(double * x, double * y, int points);
      void Polygon(double * x, double * y, int points);

      void DrawCircle(double x, double y, double r);
      void FillCircle(double x, double y, double r);
      void Circle(double x, double y, double r);

      void   SetFont(PDFFonts font, bool bold = false, bool italic = false);
      void   SetFontSize(double point_size);
      void   SetFontWidth(double relative_width);
      void   SetFontOrientation(double degrees);

      void   WriteText(double x, double y, const char * string);
      double TextWidth(const char * string);
      double TextHeight(const char * string);
      double TextExtent(const char * string);

      void   SetClipRectangle(double x0, double y0, double x1, double y1);
      void   ClearClipRectangle();

      void SetLineWidth(double width_in_points);
      void SetLineCap  (PDFLineCap lineCap);
      void SetLineStyle(PDFLineStyle lineStyle);

      void PathMoveTo(double x, double y);
      void PathLineTo(double x, double y);
      void PathBezier(double x1, double y1, double x2, double y2, double x3, double y3);
      void PathRectangle(double x1, double y1, double x2, double y2);
      void PathClose();
      void PathStroke();
      void PathFill();
      void PathStrokeAndFill();

      void SelectTextMode();
      void SelectDrawMode();

      int  GetPageNumber()   { return pages.Length(); }

      void TransformCoordinates(double xdelta, double xscale, double ydelta, double yscale);
      void RestoreCoordinates();

   private:
      IntArray pages;
      IntArray streams;
      PDF &    pdf;

      IntArray mediaBox;
      IntArray defaultBox;

      int pageRotation;
      int defaultRotation;

      bool textMode;

      int    fontId, lastFontId;
      double fontSize, lastFontSize;
      double fontWidth, lastFontWidth;
      double fontOrientation, lastFontOrientation;
      double textSin, textCos;
      double textX, textY;

      double lineRed, lineGreen, lineBlue;
      double fillRed, fillGreen, fillBlue;
      double lineCyan, lineMagenta, lineYellow, lineBlack;
      double fillCyan, fillMagenta, fillYellow, fillBlack;
      double lineGray, fillGray;
      double lastX, lastY;

      bool clip;
   };

#endif


 
 
 
