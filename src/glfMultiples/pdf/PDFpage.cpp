////////////////////////////////////////////////////////////////////// 
// pdf/PDFpage.cpp 
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
 
#include "PDFpage.h"
#include "PDF.h"

#include <math.h>

PDFPage::PDFPage(PDF & parent) : pdf(parent)
   {
   // Default to letter sized paper
   mediaBox.Dimension(4);
   mediaBox[0] = 0;
   mediaBox[1] = 0;
   mediaBox[2] = 612;
   mediaBox[3] = 792;

   defaultRotation = pageRotation = 0;

   fontId = 0;
   fontSize = 12;
   fontWidth = 1.0;

   hTextAlignment = taLeft;
   vTextAlignment = taAbove;
   }

void PDFPage::SetSize(PDFPageSize size)
   {
   switch (size)
      {
      case psLetter :
         SetSize(612, 792);
         break;
      case psLetterR :
         SetSize(792, 612);
         break;
      case psA4 :
         SetSize(595, 841);
         break;
      case psA4R :
         SetSize(841, 595);
         break;
      }
   }

void PDFPage::SetSize(int width, int height)
   {
   mediaBox[0] = 0;
   mediaBox[1] = 0;
   mediaBox[2] = width;
   mediaBox[3] = height;
   }

void PDFPage::OpenPage()
   {
   if (streams.Length())
      ClosePage();

   streams.Clear();
   streams.Push(pdf.OpenStream());

   lineRed = lineGreen = lineBlue = -1.0;
   fillRed = fillGreen = fillBlue = -1.0;
   lineCyan = lineMagenta = lineYellow = -1.0;
   fillCyan = fillMagenta = fillYellow = -1.0;
   lineGray = fillGray = -1.0;

   fontOrientation = 0.0 ;
   textCos = 1.0;
   textSin = 0.0;

   textMode = false;
   clip = false;
   }

void PDFPage::ClosePage()
   {
   if (streams.Length() == 0)
      return;

   SelectDrawMode();
   ClearClipRectangle();

   if (streams.Length())
      pdf.CloseStream();

   if (pages.Length() == 0)
      {
      defaultBox = mediaBox;
      defaultRotation = pageRotation;
      tree_index = pdf.GetObject();
      }

   int object = pdf.GetObject();
   pdf.OpenObject(object);
   pdf.OpenDictionary();

   pdf.WriteName("Type", "Page");
   pdf.WriteReference("Parent", tree_index);
   pdf.WriteReferenceArray("Contents", streams);

   pdf.WriteName("Resources");
   pdf.OpenDictionary();
   pdf.WriteName("ProcSet");
   pdf.OpenArray();
   pdf.WriteName("PDF");
   pdf.WriteName("Text");
   pdf.CloseArray();
   pdf.font.WriteResources();
   pdf.CloseDictionary();
   pdf.LineBreak();

   if (defaultBox != mediaBox)
      pdf.WriteArray("MediaBox", mediaBox);

   if (pageRotation != defaultRotation)
      pdf.WriteInteger("Rotate", pageRotation);

   pdf.CloseDictionary();
   pdf.CloseObject();

   pages.Push(object);
   streams.Clear();
   }

void PDFPage::WritePageTree()
   {
   pdf.OpenObject(tree_index);
   pdf.OpenDictionary();

   pdf.WriteName("Type", "Pages");
   pdf.WriteReferenceArray("Kids", pages);
   pdf.WriteInteger("Count", pages.Length());
   pdf.WriteArray("MediaBox", defaultBox);
   pdf.WriteInteger("Rotate", defaultRotation);

   pdf.CloseDictionary();
   pdf.CloseObject();
   }

/*
void PDFPage::AlternateZero()
   {
   if (alternate)
      return;

   pdf.AppendToStream("1 0 0 -1 0 %d cm\n", mediaBox[3]);
   SetFontOrientation(180);

   alternate = true;
   }
*/

void PDFPage::SetLineWidth(double w)
   {
   pdf.AppendToStream("%.2f w\n", w);
   }

void PDFPage::SetLineCap(PDFLineCap c)
   {
   pdf.AppendToStream("%d J\n", (int) c);
   }

void PDFPage::SetLineStyle(PDFLineStyle style)
   {
   switch (style)
      {
      case lsSolid :
         pdf.AppendToStream("[] 0 d\n");
         break;
      case lsDotted :
         pdf.AppendToStream("[0.3 3] 0 d\n");
         break;
      case lsDashed :
         pdf.AppendToStream("[3] 0 d\n");
         break;
      }
   }

void PDFPage::PathMoveTo(double x, double y)
   {
   SelectDrawMode();
   pdf.AppendToStream("%.1f %.1f m\n", x, y);
   lastX = x;
   lastY = y;
   }

void PDFPage::PathLineTo(double x, double y)
   {
   if (fabs(lastX - x) > 0.05 || fabs(lastY - y) > 0.05)
      {
      pdf.AppendToStream("%.1f %.1f l\n", x, y);
      lastX = x;
      lastY = y;
      }
   }

void PDFPage::PathBezier(double x1, double y1, double x2, double y2, double x3, double y3)
   {
   pdf.AppendToStream("%.1f %.1f %.1f %.1f %.1f %.1f c\n", x1, y1, x2, y2, x3, y3);
   }

void PDFPage::PathRectangle(double x1, double y1, double x2, double y2)
   {
   SelectDrawMode();

   if (x2 < x1) { double swap = x1; x1 = x2; x2 = swap; }
   if (y2 < y1) { double swap = y1; y1 = y2; y2 = swap; }

   pdf.AppendToStream("%.1f %.1f %.1f %.1f re\n", x1, y1, x2 - x1, y2 - y1);
   }

void PDFPage::PathClose()
   {
   pdf.AppendToStream("h\n");
   }

void PDFPage::PathStroke()
   {
   pdf.AppendToStream("S\n\n");
   }

void PDFPage::PathFill()
   {
   pdf.AppendToStream("f*\n\n");
   }

void PDFPage::PathStrokeAndFill()
   {
   pdf.AppendToStream("B*\n\n");
   }

void PDFPage::SetLineColor(double red, double green, double blue)
   {
   if (lineRed != red || lineGreen != green || lineBlue != blue)
      {
      pdf.AppendToStream("%.3f %.3f %.3f RG\n", red, green, blue);
      lineRed = red; lineBlue = blue; lineGreen = green;
      lineCyan = lineMagenta = lineYellow = lineBlack = -1.0;
      lineGray = -1.0;
      }
   }

void PDFPage::SetLineCMYK(double cyan, double magenta, double yellow, double black)
   {
   if (lineCyan != cyan || lineMagenta != magenta || lineYellow != yellow || lineBlack != black)
      {
      pdf.AppendToStream("%.3f %3f %.3f %.3f K\n", cyan, magenta, yellow, black);
      lineCyan = cyan; lineMagenta = magenta; lineYellow = yellow; lineBlack = black;
      lineRed = lineBlue = lineGreen = -1.0;
      lineGray = -1.0;
      }
   }

void PDFPage::SetLineGray(double gray)
   {
   if (lineGray != gray)
      {
      pdf.AppendToStream("%.3f G\n", gray);
      lineGray = gray;
      lineRed = lineBlue = lineGreen = -1.0;
      lineCyan = lineMagenta = lineYellow = lineBlack = -1.0;
      }
   }

void PDFPage::SetFillColor(double red, double green, double blue)
   {
   if (fillRed != red || fillGreen != green || fillBlue != blue)
      {
      pdf.AppendToStream("%.3f %.3f %.3f rg\n", red, green, blue);
      fillRed = red; fillGreen = green; fillBlue = blue;
      fillCyan = fillMagenta = fillYellow = fillBlack = -1.0;
      fillGray = -1.0;
      }
   }

void PDFPage::SetFillCMYK(double cyan, double magenta, double yellow, double black)
   {
   if (fillCyan != cyan || fillMagenta != magenta || fillYellow != yellow || fillBlack != black)
      {
      pdf.AppendToStream("%.3f %3f %.3f %.3f k\n", cyan, magenta, yellow, black);
      fillCyan = cyan; fillMagenta = magenta; fillYellow = yellow; fillBlack = black;
      fillRed = fillBlue = fillGreen = -1.0;
      fillGray = -1.0;
      }
   }

void PDFPage::SetFillGray(double gray)
   {
   if (fillGray != gray)
      {
      pdf.AppendToStream("%.3f g\n", gray);
      fillGray = gray;
      fillRed = fillBlue = fillGreen = -1.0;
      fillCyan = fillMagenta = fillYellow = fillBlack = -1.0;
      }
   }

void PDFPage::DrawLine(double x0, double y0, double x1, double y1)
   {
   PathMoveTo(x0, y0);
   PathLineTo(x1, y1);
   PathStroke();
   }

void PDFPage::DrawPolygon(double * x, double * y, int points)
   {
   PathMoveTo(x[0], y[0]);
   for (int i = 1; i < points; i++)
      PathLineTo(x[i], y[i]);
   PathClose();
   PathStroke();
   }

void PDFPage::FillPolygon(double * x, double * y, int points)
   {
   PathMoveTo(x[0], y[0]);
   for (int i = 1; i < points; i++)
      PathLineTo(x[i], y[i]);
   PathClose();
   PathFill();
   }

void PDFPage::Polygon(double * x, double * y, int points)
   {
   PathMoveTo(x[0], y[0]);
   for (int i = 1; i < points; i++)
      PathLineTo(x[i], y[i]);
   PathClose();
   PathStrokeAndFill();
   }

void PDFPage::DrawRectangle(double x0, double y0, double x1, double y1)
   {
   PathRectangle(x0, y0, x1, y1);
   PathStroke();
   }

void PDFPage::FillRectangle(double x0, double y0, double x1, double y1)
   {
   PathRectangle(x0, y0, x1, y1);
   PathFill();
   }

void PDFPage::Rectangle(double x0, double y0, double x1, double y1)
   {
   PathRectangle(x0, y0, x1, y1);
   PathStrokeAndFill();
   }

void PDFPage::DrawCircle(double x, double y, double r)
   {
   double k = 0.5522847498 * r;

   PathMoveTo(x, y + r);
   PathBezier(x + k, y + r, x + r, y + k, x + r, y);
   PathBezier(x + r, y - k, x + k, y - r, x, y - r);
   PathBezier(x - k, y - r, x - r, y - k, x - r, y);
   PathBezier(x - r, y + k, x - k, y + r, x, y + r);
   PathClose();
   PathStroke();
   }

void PDFPage::FillCircle(double x, double y, double r)
   {
   double k = 0.5522847498 * r;

   PathMoveTo(x, y + r);
   PathBezier(x + k, y + r, x + r, y + k, x + r, y);
   PathBezier(x + r, y - k, x + k, y - r, x, y - r);
   PathBezier(x - k, y - r, x - r, y - k, x - r, y);
   PathBezier(x - r, y + k, x - k, y + r, x, y + r);
   PathClose();
   PathFill();
   }

void PDFPage::Circle(double x, double y, double r)
   {
   double k = 0.5522847498 * r;

   PathMoveTo(x, y + r);
   PathBezier(x + k, y + r, x + r, y + k, x + r, y);
   PathBezier(x + r, y - k, x + k, y - r, x, y - r);
   PathBezier(x - k, y - r, x - r, y - k, x - r, y);
   PathBezier(x - r, y + k, x - k, y + r, x, y + r);
   PathClose();
   PathStrokeAndFill();
   }

void PDFPage::SelectTextMode()
   {
   if (textMode)
      return;

   pdf.AppendToStream("BT\n");

/*
   if (alternate)
      pdf.AppendToStream("1 0 0 -1 0 0 Tm\n");
*/
   textMode = true;

   lastFontId = -1;
   lastFontSize = -1;
   lastFontWidth = 1.0;

   lastFontOrientation = 0;

   textX = textY = 0;
   }

void PDFPage::SelectDrawMode()
   {
   if (!textMode)
      return;

   pdf.AppendToStream("ET\n");

   textMode = false;
   }

double PDFPage::TextWidth(const char * string)
   {
   return fabs(TextExtent(string) * textCos - fontSize * textSin);
   }

double PDFPage::TextHeight(const char * string)
   {
   return fabs(fontSize * textCos + TextExtent(string) * textSin);
   }

double PDFPage::TextExtent(const char * string)
   {
   return pdf.font.TextWidth(fontId, string) * fontSize * 0.001 * fontWidth;
   }

void PDFPage::WriteText(double x, double y, const char * text)
   {
   SelectTextMode();

   if (fontId != lastFontId || fontSize != lastFontSize)
      pdf.font.SelectFont(lastFontId = fontId, lastFontSize = fontSize);

   if (fontWidth != lastFontWidth)
      pdf.AppendToStream("%.1f Tz\n", (lastFontWidth = fontWidth) * 100.);

   if (fontOrientation != lastFontOrientation)
      {
      pdf.AppendToStream("%.3f %.3f %.3f %.3f 0 0 Tm\n",
                         textCos, textSin, -textSin, textCos);

      textX = textY = 0;

      lastFontOrientation = fontOrientation;
      }

   // Calculate bounding rectangle for text
   double extent = TextExtent(text);

   double box_x = textCos * extent - textSin * fontSize;
   double box_y = textCos * fontSize + textSin * extent;

   // Adjust horizontal placement
   switch (hTextAlignment)
      {
      case taLeft :
         if (box_x < 0) x -= box_x;
         break;
      case taRight :
         if (box_x > 0) x -= box_x;
         break;
      case taCenter :
         x -= box_x * 0.5;
         break;
      }

   // Adjust vertical placement
   switch (vTextAlignment)
      {
      case taAbove :
         if (box_y < 0) y -= box_y;
         break;
      case taBelow :
         if (box_y > 0) y -= box_y;
         break;
      case taMiddle :
         y -= box_y * 0.5;
         break;
      }

   // Move the text cursor to appropriate position, taking transformation
   // matrix into account ...
   double deltaX = x - textX;
   double deltaY = y - textY;

   double textModeX = deltaX * textCos + deltaY * textSin;
   double textModeY = deltaY * textCos - deltaX * textSin;

   pdf.AppendToStream("%.2f %.2f Td\n", textModeX, textModeY);
   pdf.WriteString(text);
   pdf.AppendToStream(" Tj\n");

   textX = x;
   textY = y;
   }

void PDFPage::SetFont(PDFFonts font, bool bold, bool italic)
   {
   fontId = pdf.font.GetFontID(font, bold, italic);
   }

void PDFPage::SetFontSize(double pointSize)
   {
   fontSize = pointSize;
   }

void PDFPage::SetFontWidth(double width)
   {
   fontWidth = width;
   }

void PDFPage::SetFontOrientation(double orientation)
   {
   double rad = orientation * 0.0174532925199;

   textCos = cos(rad);
   textSin = sin(rad);

   fontOrientation = orientation;
   }

void PDFPage::SetClipRectangle(double x0, double y0, double x1, double y1)
   {
   SelectDrawMode();
   ClearClipRectangle();

   pdf.AppendToStream("q\n");
   PathRectangle(x0, y0, x1, y1);
   pdf.AppendToStream("W n\n\n");

   clip = true;
   }

void PDFPage::ClearClipRectangle()
   {
   if (clip)
      {
      SelectDrawMode();
      pdf.AppendToStream("Q\n\n");
      clip = false;
      }
   }

double PDFPage::GetHeight()
   {
   return fabs((double) (mediaBox[3] - mediaBox[1]));
   }

double PDFPage::GetWidth()
   {
   return fabs((double) (mediaBox[2] - mediaBox[0]));
   }


void PDFPage::TransformCoordinates(double xdelta, double xscale, double ydelta, double yscale)
   {
   SelectDrawMode();
   pdf.AppendToStream("q\n");
   pdf.AppendToStream("%g 0 0 %g %g %g cm\n", xscale, yscale, xdelta, ydelta);
   }

void PDFPage::RestoreCoordinates()
   {
   SelectDrawMode();
   pdf.AppendToStream("Q\n");
   }






 
