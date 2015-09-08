#include "JevpPlot.h"
#include <TLegend.h>
#include <TCanvas.h>
#include <TStyle.h>

PlotHisto::PlotHisto()
{
  histo = NULL;
  legendText = NULL;
  legendArgs = NULL;
}

PlotHisto::PlotHisto(PlotHisto &x)
{
  if(x.histo) {
    histo = (TH1 *)x.histo->Clone();
  }
  else 
    histo = NULL;

  if(x.legendText) {
    legendText = new char[sizeof(x.legendText) + 1];
    strcpy(legendText, x.legendText);
  }
  else 
    legendText = NULL;

  if(x.legendArgs) {
    legendArgs = new char[sizeof(x.legendArgs) + 1];
    strcpy(legendArgs, x.legendArgs);
  }
  else 
    legendArgs = NULL;
}


void PlotHisto::setLegText(char *text)
{
  if(legendText) delete legendText;
  legendText = new char[strlen(text)+1];
  strcpy(legendText, text);
}

void PlotHisto::setLegArgs(char *text)
{
  if(legendArgs) delete legendText;
  legendArgs = new char[strlen(text)+1];
  strcpy(legendArgs, text);
}

// PlotHisto::PlotHisto(TH1 *hist, char *legText, char *legArgs) {
//   histo = hist;

//   if(!legText) legendText = NULL;
//   else {
//     legendText = new char[strlen(legText)+1];
//     strcpy(legendText, legText);
//   }
  
//   if(!legArgs) legendArgs = NULL;
//   else {
//     legendArgs = new char[strlen(legArgs)+1];
//     strcpy(legendArgs, legArgs);
//   }
// }

PlotHisto::~PlotHisto()
{
  if(histo) delete histo;
  if(legendText) delete legendText;
  if(legendArgs) delete legendArgs;
}


PlotTF1::PlotTF1()
{
	histo = NULL;
	legendText = NULL;
	legendArgs = NULL;
}

PlotTF1::PlotTF1(PlotTF1 &x)
{
	if(x.histo) {
		histo = (TF1 *)x.histo->Clone();
	}
	else 
		histo = NULL;

	if(x.legendText) {
		legendText = new char[sizeof(x.legendText) + 1];
		strcpy(legendText, x.legendText);
	}
	else 
		legendText = NULL;

	if(x.legendArgs) {
		legendArgs = new char[sizeof(x.legendArgs) + 1];
		strcpy(legendArgs, x.legendArgs);
	}
	else 
		legendArgs = NULL;
}


void PlotTF1::setLegText(char *text)
{
	if(legendText) delete legendText;
	legendText = new char[strlen(text)+1];
	strcpy(legendText, text);
}

void PlotTF1::setLegArgs(char *text)
{
	if(legendArgs) delete legendText;
	legendArgs = new char[strlen(text)+1];
	strcpy(legendArgs, text);
}

// PlotTF1::PlotTF1(TH1 *hist, char *legText, char *legArgs) {
//   histo = hist;

//   if(!legText) legendText = NULL;
//   else {
//     legendText = new char[strlen(legText)+1];
//     strcpy(legendText, legText);
//   }

//   if(!legArgs) legendArgs = NULL;
//   else {
//     legendArgs = new char[strlen(legArgs)+1];
//     strcpy(legendArgs, legArgs);
//   }
// }

PlotTF1::~PlotTF1()
{
	if(histo) delete histo;
	if(legendText) delete legendText;
	if(legendArgs) delete legendArgs;
}


void JevpPlot::setDrawOpts(char *opts)
{
  if(drawopts) delete drawopts;
  drawopts = new char[strlen(opts)+1];
  strcpy(drawopts, opts);
}

JevpPlot::JevpPlot() { 
  logx=0;
  logy=0;
  nhistos=0;
  legendx1 = -1;
  legendx2 = -1;
  legendy1 = -1;
  legendy2 = -1;
  optstat = 1;
  drawopts = NULL;
  xaxis_label = NULL;
  yaxis_label = NULL;
  legend = NULL;

  run = 0;
  refid = 0;
  refcomment = NULL;
}

JevpPlot::JevpPlot(JevpPlot &x)
{
  if(x.legend) {
    legend = new TLegend(*(x.legend));
  }
  else legend = NULL;
  
  run = x.run;
  refid = x.refid;
  if(x.refcomment) {
    refcomment = new char[strlen(x.refcomment)+1];
    strcpy(refcomment, x.refcomment);
  }
  else
    refcomment = NULL;

  optstat = x.optstat;
  logx = x.logx;
  logy = x.logy;
  legendx1 = x.legendx1;
  legendx2 = x.legendx2;
  legendy1 = x.legendy1;
  legendy2 = x.legendy2;
  if(x.xaxis_label) {
    xaxis_label = new char[strlen(x.xaxis_label) + 1];
    strcpy(xaxis_label, x.xaxis_label);
  }
  else
    xaxis_label = NULL;
  if(x.yaxis_label) {
    yaxis_label = new char[strlen(x.yaxis_label) + 1];
    strcpy(yaxis_label, x.yaxis_label);
  }
  else
    yaxis_label = NULL;
  
  PlotHisto *curr = (PlotHisto *)x.histos.First();
  while(curr) {
    PlotHisto *nph = new PlotHisto(*curr);
    histos.Add(nph);
    curr = (PlotHisto *)x.histos.After(curr);
  }

  nhistos = x.nhistos;
  if(x.drawopts) {
    drawopts = new char[strlen(x.drawopts) + 1];
    strcpy(drawopts, x.drawopts);
  }
  else
    drawopts = NULL;
}


void JevpPlot::setRefComment(char *text)
{
  if(refcomment) delete refcomment;
  refcomment = new char[strlen(text)+1];
  strcpy(refcomment, text);
}

// Carefull!  This deletes the memory for the histograms...
// A strange usage, but neccessary because the JevpPlot is
// transferred via ethernet to processes that know nothing
// about the original PlotHisto memory allocations...
//
// a.  histos added must be dynamically created using "new"
// b.  histos memory is owned by JevpPlot once added
JevpPlot::~JevpPlot() {  


  //printf("Remove the histos.... ~JevpPlot\n");
  PlotHisto *curr = (PlotHisto *)histos.First();

  printf("DESTRUCTOR FOR PLOT:  %s %d\n",(curr != NULL) ? curr->GetName() : "noname", refid);

  PlotHisto *next;
  while(curr) {
    next = (PlotHisto *)histos.After(curr);
    histos.Remove(curr);  // remove from list
    delete curr;          // delete

    curr = next;
  }

  if(drawopts) {
    //printf("Drawopts...\n");
    delete drawopts;
  }
  
  if(legend) {
    printf("Delete legend\n");
    delete legend;
    printf("Done delete legend\n");
  }
  
  if(refcomment) {
    delete refcomment;
  }

  printf("Done with jevpplot destructor...\n");
}


void JevpPlot::addHisto(PlotHisto *hist) {
  histos.Add(hist);
  nhistos++;
}

void JevpPlot::addHisto(PlotTF1 *hist){
histos.Add(hist);
nhistos++;
}


PlotHisto *JevpPlot::getHisto(int i) {
  int idx=0;
  
  PlotHisto *curr = (PlotHisto *)histos.First();

  while(curr) {
    if(i == idx) return curr;
    idx++;
    curr = (PlotHisto *)histos.After(curr);
  }
  
  return NULL;
}

void JevpPlot::removeHisto(int i)
{
  int idx = 0;
  PlotHisto *curr = (PlotHisto *)histos.First();
  
  while(curr) {
    if(idx == i) {
      histos.Remove(curr);
      nhistos--;
      return;
    }

    curr = (PlotHisto *)histos.After(curr);
  }
}


int JevpPlot::nHistos()
{
  return nhistos;
}

char *JevpPlot::GetPlotName()
{
  PlotHisto *curr = (PlotHisto *)histos.First();

  if(!curr) return NULL;
  
  if(!curr->histo) return NULL;

  return (char *)curr->histo->GetName();
}


void JevpPlot::draw()
{
	// Check for legends...
	//  PlotHisto *curr = (PlotHisto *)histos.First();

	TObject* curr = histos.First();

	if(legend != NULL) {
		printf("draw delete legend\n");
		delete legend;
		printf("draw delete legend done\n");
		legend = NULL;
	}

	if(legendx1 >= 0) {
		legend = new TLegend(legendx1,legendy1,legendx2,legendy2);
	}

	gStyle->SetOptStat(optstat);
	gPad->SetLogx(logx);
	gPad->SetLogy(logy);

	char *same = NULL;

	while(curr){
		if (curr->InheritsFrom( PlotHisto::Class())) {

			PlotHisto* currHisto = (PlotHisto *)curr;

			if(logy) {
				if(currHisto->histo->GetMaximum() == 0.0) currHisto->histo->SetMaximum(10.0);
				if(currHisto->histo->GetMinimum() == 0.0) currHisto->histo->SetMinimum(1.0);
			}
			else {
				if(currHisto->histo->GetMaximum() == 0.0) currHisto->histo->SetMaximum(1.0);
			}

			if(legend) {
				char *text = (char *)((currHisto->legendText) ? currHisto->legendText : "no text");
				char *args = (char *)((currHisto->legendArgs) ? currHisto->legendArgs : "");

				legend->AddEntry(currHisto->histo, text, args); 
			}

			char opts[256];
			if(drawopts)
			strcpy(opts, drawopts);
			else opts[0] = '\0';

			if(same) strcat(opts,"SAME");

			printf("opts---%s\n",opts);
			currHisto->histo->Draw(opts);
			same = "SAME";
//			printf("same = %s",same) ;


		} 
		else if (curr->InheritsFrom( PlotTF1::Class())){
			PlotTF1* currTF1 = (PlotTF1 *)curr;

//			char opts[256];
//			if(drawopts)
//				strcpy(opts, drawopts);
//			else opts[0] = '\0';

//			printf("opts---%s\n",opts);
			//currTF1->histo->Draw(opts);
//			printf("testing draw \n") ;
			printf("opts---%s\n","SAME");
			currTF1->histo->SetLineWidth(0.5) ;
			currTF1->histo->Draw("SAME");

			//do whatever you want to do with TF1
		}
		curr = (TObject *)histos.After(curr);
	}
	if(legend)
		legend->Draw();
}

