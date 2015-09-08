#include <stdio.h>
#include <stdlib.h>

#include "JevpPlotSet.h"
#include "DAQ_READER/daqReader.h"

#include <TH1I.h>
#include <TH2F.h>

#include <math.h>

class ExampleHist : public JevpPlotSet {
public:
  JevpPlot *p1;
  JevpPlot *p2;
  JevpPlot *p3;

  void initialize(int argc, char *argv[]) {
    p1 = new JevpPlot();
    p1->logy=1;

    printf("a\n");
    TH1 *h = new TH1I("h77_zdc_time_west", "bbb",32,0,31);
    PlotHisto *ph = new PlotHisto();
    ph->histo = h;
    h->SetFillColor(5);
    printf("b\n");
    p1->addHisto(ph);
    printf("c\n");
    addPlot(p1);
    // access h1 as:  TH1 *h = getPlot(0)->getHisto(0)->histo;
    
    printf("n\n");
    
    p2 = new JevpPlot();
    p2->setLegend(.5,.5,.7,.7);
    p2->optstat = 0;

    h = new TH1F("h76_zdc_time_east","aaa",21,0,21);
    h->SetFillColor(1);
    h->SetLineColor(1);
    ph = new PlotHisto();
    ph->histo = h;
    ph->setLegText("legend1");
    ph->setLegArgs("F");
    p2->addHisto(ph);

    h = new TH1F("blah","bc1-a",21,0,21);
    h->SetFillColor(3);
    h->SetLineColor(3);
    ph = new PlotHisto();  
    ph->histo = h;
    ph->setLegText("legend2");
    ph->setLegArgs("F");
    p2->addHisto(ph);
    
    addPlot(p2);

    p3 = new JevpPlot();
    p3->optstat = 0;
    p3->setDrawOpts("colz");
    h = new TH2F("h78_zdc_timediff_east_west","ccc",120,0,240,100,0,15);
    ph = new PlotHisto();
    ph->histo = h;
    p3->addHisto(ph);
    addPlot(p3);
  };
  
  void startrun(daqReader *rdr) {
    getPlotByIndex(0)->getHisto(0)->histo->Reset();
    getPlotByIndex(1)->getHisto(0)->histo->Reset();
    getPlotByIndex(1)->getHisto(1)->histo->Reset();
    printf("Starting run #%d\n",rdr->run);
    
  };

  void stoprun(daqReader *rdr) {
    printf("Stopping run #%d\n",rdr->run);
  };
  
  void event(daqReader *rdr) {
    static int tm=0;
    
    //sleep(1);

    u_int trg = rdr->daqbits;
    printf("Event %d/%d: 0x%x\n",rdr->event_number,rdr->seq,trg);
    for(int i=0;i<32;i++) {
      if(trg & (1<<i)) {
	((TH1I *)getPlotByIndex(0)->getHisto(0)->histo)->Fill(i);
      }
    }
    ((TH1F *)getPlotByIndex(1)->getHisto(0)->histo)->Fill(10);
    ((TH1F *)getPlotByIndex(1)->getHisto(1)->histo)->Fill(15);
    
    
    if(tm ==0) tm = rdr->evt_time;

    printf("sz=%d %d\n",rdr->event_size,rdr->evt_time-tm);
    ((TH2F *)getPlotByIndex(2)->getHisto(0)->histo)->Fill(rdr->evt_time - tm,log(rdr->event_size));
  };

  int selectEvent(daqReader *rdr) {
    return 1;
  };

  int selectRun(daqReader *rdr) {
    return 1;
  };  
};


int main(int argc, char *argv[])
{
  ExampleHist me;
  
  me.Main(argc, argv);
}
