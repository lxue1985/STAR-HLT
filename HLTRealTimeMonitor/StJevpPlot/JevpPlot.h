#ifndef _JEVPPLOT_H_
#define _JEVPPLOT_H_

#include <TROOT.h>
#include <TH1.h>
#include <TF1.h>
#include <TObject.h>
#include <TLegend.h>

// This class contains everything needed to draw a 
// single plot, whether it be a single histogram
// or a complex plot made from multiple histograms with 
// multiple colors/legends/whatever...

class PlotHisto : public TObject {
 public:
  TH1 *histo;

  char *legendText;
  char *legendArgs;
  
  PlotHisto();
  PlotHisto(PlotHisto &x);

  PlotHisto(const PlotHisto &x) {
    printf("Copy constructor\n");
    histo = new TH1(*(x.histo));
    if(x.legendText) setLegText(x.legendText);
    if(x.legendArgs) setLegArgs(x.legendArgs);
  }

  // PlotHisto(TH1 *hist, char *legText=NULL, char *legArgs=NULL);
  void setLegText(char *text);
  void setLegArgs(char *text);

  virtual ~PlotHisto();

  ClassDef(PlotHisto, 1);
};

class PlotTF1 : public TObject {
	public:
		TF1 *histo;

		char *legendText;
		char *legendArgs;

		PlotTF1();
		PlotTF1(PlotTF1 &x);

		PlotTF1(const PlotTF1 &x) {
			printf("Copy constructor\n");
			histo = new TF1(*(x.histo));
			if(x.legendText) setLegText(x.legendText);
			if(x.legendArgs) setLegArgs(x.legendArgs);
		}

		// PlotTF1(TH1 *hist, char *legText=NULL, char *legArgs=NULL);
		void setLegText(char *text);
		void setLegArgs(char *text);

		virtual ~PlotTF1();

		ClassDef(PlotTF1, 1);
};



class JevpPlot : public TObject {
 public:
  // int run;
  TLegend *legend;
 
  int run;       
  int refid;         // reference name, if any
  char *refcomment;    // reference comment, if any

  int optstat;
  int logx;
  int logy;
  double legendx1;
  double legendy1;
  double legendx2;
  double legendy2;

  void setRefComment(char *comment);
  void addHisto(PlotHisto *hist); 
  void addHisto(PlotTF1 *hist);

  void removeHisto(int i);
  PlotHisto *getHisto(int i);
  int nHistos();
  char *GetPlotName();
  void draw();
  void setLegend(double x1,double y1, double x2, double y2) {
    legendx1=x1;legendx2=x2;legendy1=y1;legendy2=y2;
  };
  void setDrawOpts(char *opts);
  JevpPlot();
  JevpPlot(JevpPlot &x);
  ~JevpPlot();
  
 private:
  char *xaxis_label;
  char *yaxis_label;

  TList histos;
  int nhistos;
  char *drawopts;
  ClassDef(JevpPlot, 1);
};


#endif
