using namespace RooFit;

EColor csig = kRed;
EColor csi2 = kViolet;
EColor cbkg = kGreen;
EColor cmis = kOrange;

#define SETLOGVAL(X,Y) bool X = ( Y.find(#X) == Y.end() ? false : true );

/*
bool fitijp = false;
bool fitbjp = false;
bool fitphi = false;
bool fitip2 = false;
bool fitbp2 = false;
bool fitups = false;
bool fitbui = false;
bool fitbud = false;
bool fitbsd = false;
bool fitbdd = false;
bool fitb0i = false;
bool fitlbi = false;
bool fitlbd = false;
bool fitbci = false;
bool fitbcd = false;
bool fitx3i = false;
bool fitx3d = false;

bool fitA = false;
bool fitB = false;
bool fitC = false;
*/

set<string> parseLabels( const string& labelString ) {
  set<string> labelSet;
//  labelSet.clear();
  const char* labelPtr = labelString.c_str();
  unsigned long labelLength = labelString.length();
  while ( labelLength ) {
    while ( *labelPtr == ' ' ) ++labelPtr;
    string labelNext( labelPtr );
    cout << "parse ---" << labelNext << "---" << endl;
    if ( labelNext == "" ) break;
    unsigned long currentLength = labelLength = labelNext.find( ":" );
    string labelCurr ;
    if ( currentLength != string::npos ) {
      labelCurr = labelNext.substr( 0, currentLength );
      labelPtr += ++currentLength;
    }
    else {
      labelCurr = labelNext;
      labelLength = 0;
    }
    currentLength = labelCurr.find( " " );
    if ( !currentLength ) break;
    if ( currentLength != string::npos )
         labelCurr = labelCurr.substr( 0, currentLength );
    labelSet.insert( labelCurr );
  }
  return labelSet;
}

double errFunct( double x, double threshold, double slope ) {
  float mmin = 5.0;
  float mmax = 6.0;
/*
class ErrF {
 private:
  float mmin;
  float mmax;
 public:
  ErrF( float massMin, float massMax ): mmin( massMin ), mmax( massMax ) {}
  ~ErrF() {}
  float operator()( double x, double threshold, double slope ) const {
*/
    static TF1* tferf = new TF1( "tferf", "1-erf((x-[0])/[1])", mmin, mmax );
  tferf->SetParameter( 0, threshold );
  tferf->SetParameter( 1, slope );
  return ( *tferf )( x );
}
//};

void print( const RooRealVar& x, ofstream& os, float s = 1.0 ) {
  cout << x.GetName();
  os   << '$'
       << x.GetName();
  if ( s != 1.0 ) {
    cout << "*" << s;
    os   << "*" << s;
  }
  cout << " = " << scientific << setw( 20 ) << setprecision( 12 )
       << x.getVal()   * s << " +/- "
       << x.getError() * s << endl;
  os   << " = " << scientific << setw( 20 ) << setprecision( 12 )
       << x.getVal  () * s << " \\pm "
       << x.getError() * s << '~'
       << x.getUnit ()     << '$' << endl;
  os   << "%" << scientific << setw( 20 ) << setprecision( 12 )
       << x.getVal  () * s << " +- "
       << x.getError() * s << '~'
       << x.getUnit ()     << endl;
  return;
}

float linTrans( float min, float max, float fac ) {
  return ( min * ( 1.0 - fac ) ) +
         ( max *         fac   );
}

float linTrans( RooPlot* frame, float fac ) {
  return linTrans( frame->GetMinimum(), frame->GetMaximum(), fac );
}

float logTrans( RooPlot* frame, float fac ) {
  return pow( frame->GetMinimum(), 1.0 - fac ) *
         pow( frame->GetMaximum(),       fac );
}

void printSigma( RooPlot* frame, float sigma, float xmin, float xmax,
                 float xfac = 0.65, float yfac = 0.85 ) {
  TLatex t;
  string s = "#font[62]{#sigma = ";
  s += to_string( lround( sigma * 1000 ) );
  s += " MeV}";
  t.SetTextSize( frame->GetYaxis()->GetLabelSize() );
//  t.DrawLatex( ( xmin * ( 1.0 - xfac ) ) + ( xmax * xfac ),
  t.DrawLatex( linTrans( xmin, xmax, xfac ),
               frame->GetMaximum() * yfac, s.c_str() );
}

class InfoLabel {
 private:
  float yfac;
  string label;
 public:
  enum mode { lin, log };
  InfoLabel( const string& s ): label( s ) {}
  ~InfoLabel() {}
  void print( RooPlot* frame,
              float xmin, float xmax,
              float xfac = 0.65, float yfac = 0.85, mode m = lin ) const;
/*
              float xfac = 0.65, float yfac = 0.85, mode m = lin ) const {
    TLatex t;
    t.SetTextSize( frame->GetYaxis()->GetLabelSize() );
    switch ( m ) {
    case lin:
//      cout << "PRINT INFO " << m << ' ' << xfac << ' ' << linTrans( xmin, xmax, xfac ) << ' ' << yfac << ' ' << logTrans( frame, yfac ) << endl;
//      t.DrawLatex( ( xmin * ( 1.0 - xfac ) ) + ( xmax * xfac ),
//                   frame->GetMaximum() * yfac, label.c_str() );
      t.DrawLatex( linTrans( xmin, xmax, xfac ),
                   linTrans( frame, yfac ), label.c_str() );
      break;
    case log:
//      cout << "PRINT INFO " << m << ' ' << xfac << ' ' << linTrans( xmin, xmax, xfac ) << ' ' << yfac << ' ' << logTrans( frame, yfac ) << endl;
      t.DrawLatex( linTrans( xmin, xmax, xfac ),
                   logTrans( frame, yfac ), label.c_str() );
//      t.DrawLatex( ( xmin * ( 1.0 - xfac ) ) + ( xmax * xfac ),
//                   frame->GetMaximum() * yfac, label.c_str() );
      break;
    }
    return;
  }
*/
};

float getYPos( RooPlot* frame, float yfac,
               InfoLabel::mode m = InfoLabel::lin ) {
  switch ( m ) {
  case InfoLabel::lin: return linTrans( frame, yfac );
  case InfoLabel::log: return logTrans( frame, yfac );
  }
  return 0.0;
}

void InfoLabel::print( RooPlot* frame,
                      float xmin, float xmax,
                      float xfac = 0.65, float yfac = 0.85,
                      InfoLabel::mode m = InfoLabel::lin ) const {
  TLatex t;
  t.SetTextSize( frame->GetYaxis()->GetLabelSize() );
  t.DrawLatex( linTrans( xmin, xmax, xfac ),
               getYPos( frame, yfac, m ), label.c_str() );
  return;
}

void printInfo( RooPlot* frame, const vector<const InfoLabel*>& lv,
                float xmin, float xmax,
                float xfac = 0.65, float yfac = 0.85,
                InfoLabel::mode m = InfoLabel::lin ) {
  int nl = lv.size();
  int il;
  for ( il = 0; il < nl; ++il ) lv[il]->print( frame, xmin, xmax,
                                               xfac, yfac -= 0.06, m );
}

void printLine( RooPlot* frame,
                float xmin, float xmax, float xfac, float yfac,
                const string& label,
                ELineStyle style = kSolid, EColor color = kBlue,
                InfoLabel::mode m = InfoLabel::lin ) {
  TLine l;
  l.SetLineStyle( style );
  l.SetLineColor( color );
  l.SetLineWidth( frame->GetYaxis()->GetLabelSize() * 100 );
  float ypos = getYPos( frame, yfac + 0.013, m );
  cout << "PRINT LINE " << m << ' ' << linTrans( xmin, xmax, xfac ) << ' ' << ypos << endl;
  l.DrawLine( linTrans( xmin, xmax, xfac        ), ypos, 
              linTrans( xmin, xmax, xfac + 0.07 ), ypos );
  TLatex t;
  t.SetTextSize( frame->GetYaxis()->GetLabelSize() / 1.2 );
  cout << "PRINT LINE " << frame->GetYaxis()->GetLabelSize() << endl;
  t.SetTextFont( frame->GetYaxis()->GetLabelFont() );
  t.DrawLatex( linTrans( xmin, xmax, xfac + 0.085 ),
               getYPos( frame, yfac, m ) , label.c_str() );
  return;
}
                

void printTitle( RooPlot* frame, int year, float lumi,
                 float xmin, float xmax, float xfac, float yfac,
                 InfoLabel::mode m ) {
  float ypos = getYPos( frame, yfac, m );
//  switch ( m ) {
//  case InfoLabel::lin: ypos = linTrans( frame, yfac ); break;
//  case InfoLabel::log: ypos = logTrans( frame, yfac ); break;
//  }
//  cout << "SET TITLE " << m << ' ' << xfac << ' ' << linTrans( xmin, xmax, xfac ) << ' ' << yfac << ' ' << ypos << endl;
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( frame->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  t.SetTextSize( frame->GetYaxis()->GetLabelSize() );
  t.SetTextFont( frame->GetYaxis()->GetLabelFont() );
  t.DrawLatex( linTrans( xmin, xmax, xfac ), ypos, s.c_str() );
//  t.DrawLatex( ( xmin * ( 1.0 - xfac ) ) +
//               ( xmax *         xfac   ), ypos, s.c_str() );
  xfac -= 0.35;
  t.DrawLatex( linTrans( xmin, xmax, xfac ), ypos,
//  t.DrawLatex( ( xmin * ( 1.0 - xfac ) ) +
//               ( xmax *         xfac   ), ypos,
               "#font[62]{CMS} #font[52]{Preliminary}" );
  return;
}

void printTitleLin( RooPlot* frame, int year, float lumi,
                    float xmin, float xmax,
                    float xfac = 0.55, float yfac = 1.01 ) {
  printTitle( frame, year, lumi, xmin, xmax, xfac, yfac, InfoLabel::lin );
/*
  float ypos = linTrans( frame, yfac );
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "L=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  t.SetTextSize( frame->GetYaxis()->GetLabelSize() );
  t.SetTextFont( frame->GetYaxis()->GetLabelFont() );
  t.DrawLatex( linTrans( xmin, xmax, xfac ), ypos, s.c_str() );
//  t.DrawLatex( ( xmin * ( 1.0 - xfac ) ) +
//               ( xmax *         xfac   ), ypos, s.c_str() );
  xfac -= 0.35;
  t.DrawLatex( linTrans( xmin, xmax, xfac ), ypos,
//  t.DrawLatex( ( xmin * ( 1.0 - xfac ) ) +
//               ( xmax *         xfac   ), ypos,
               "#font[62]{CMS} #font[52]{Preliminary}" );
*/
  return;
}

void printTitleLog( RooPlot* frame, int year, float lumi,
                    float xmin, float xmax,
                    float xfac = 0.55, float yfac = 1.01 ) {
  printTitle( frame, year, lumi, xmin, xmax, xfac, yfac, InfoLabel::log );
/*
//                  float xfac = 0.55, float yfac = 1.01,
//                  float topMargin = 1.2, float range = 2500.0 ) {
//  float ymax = topMargin * frame->GetMaximum();
//  frame->SetMaximum( ymax );
//  frame->SetMinimum( ymax / range );
//  float ypos = pow( frame->GetMinimum(), 1.0 - yfac ) *
//               pow( frame->GetMaximum(),       yfac );
  float ypos = logTrans( frame, yfac );
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( frame->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  t.SetTextSize( frame->GetYaxis()->GetLabelSize() );
  t.SetTextFont( frame->GetYaxis()->GetLabelFont() );
  t.DrawLatex( linTrans( xmin, xmax, xfac ), ypos, s.c_str() );
//  t.DrawLatex( ( xmin * ( 1.0 - xfac ) ) +
//               ( xmax *         xfac   ), ypos, s.c_str() );
  xfac -= 0.35;
  t.DrawLatex( linTrans( xmin, xmax, xfac ), ypos,
//  t.DrawLatex( ( xmin * ( 1.0 - xfac ) ) +
//               ( xmax *         xfac   ), ypos,
               "#font[62]{CMS} #font[52]{Preliminary}" );
*/
  return;
}

void fitCBPol2 ( const string& head,
                 const string& name,
                 const string& pset, int year,
                 const string& dset, float lumi,
                 const string& ptit,
                 const string& atit,
                 const vector<const InfoLabel*>& lv,
                 int nbin, double xmin, double xmax, float xfac, float yfac,
                 double mass, double sigm,
                 double tcut, double tpow, double fsig,
                 bool wwid = true, bool pbkg = true, bool psig = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar rmass( "rmass", "Resonance  Mass" , mass, xmin, xmax, "GeV" );
  RooRealVar width( "width", "Resonance Width" , sigm, sigm/10, sigm*2, "GeV" );
  RooRealVar sfcut( "sfcut", "Cut Distance"    , tcut, 0.0, tcut*2, "GeV" );
  RooRealVar power( "power", "Tail Power"      , tpow, 0.0, tpow*10, "." );
  RooRealVar coef1( "coef1", "Coefficicent 1"  , 0.0, -100.0, 100.0, "GeV^-1" );
  RooRealVar coef2( "coef2", "Coefficicent 1"  , 0.0,  -10.0,  10.0, "GeV^-2" );
//  RooRealVar frac1( "frac1", "signal fraction1",  0.5, 0.0, 1.0 );
//  RooRealVar frac2( "frac2", "signal fraction2", fsig, 0.0, 1.0 );
  // create the fit model components: CB and polynomial PDFs
  RooCBShape crybf( "crybf", "Crystal Ball Lineshape", m, rmass, width,
                                                          sfcut, power );
  RooPolynomial polyn( "polyn", "Polynomial Distribution", m, RooArgList( coef1, coef2 ) );
  // combine them using addition
//  RooAddPdf cbsum( "cbsum", "Signal", RooArgList( cryb1, cryb2 ), frac1 );
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, m );
  int nevr = data->numEntries();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs / 10.0, nevs * 10.0 );
  RooExtendPdf exsig( "exsig", "esig", crybf, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", polyn, nbkg );
  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbkg ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( cbsum, polyn ), frac2 );
  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data, Extended( kTRUE ) );
//  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
/*
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "width: " << width.getVal() << " +/- " << width.getError() << endl;
  cout << "sfcut: " << sfcut.getVal() << " +/- " << sfcut.getError() << endl;
  cout << "power: " << power.getVal() << " +/- " << power.getError() << endl;
  cout << "coef1: " << coef1.getVal() << " +/- " << coef1.getError() << endl;
  cout << "coef2: " << coef2.getVal() << " +/- " << coef2.getError() << endl;
//  cout << "frac1: " << frac1.getVal() << " +/- " << frac1.getError() << endl;
//  cout << "frac2: " << frac2.getVal() << " +/- " << frac2.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
*/
//  float f1 = frac1.getVal();
//  float sigma = sqrt( ( pow( widt1.getVal(), 2 ) *         f1   ) +
//                      ( pow( widt2.getVal(), 2 ) * ( 1.0 - f1 ) ) );
  float sigma = width.getVal();
  ofstream os( ( nfit + ".tex" ).c_str() );
  print( rmass, os );
  print( width, os );
  print( sfcut, os );
  print( power, os );
  print( coef1, os );
  print( coef2, os );
//  print( frac1, os );
  print( nsig , os );
  print( nbkg , os );
  print( nsig , os, 1.0 / lumi );
  cout << " sigma = " << sigma << endl;
  os << "$\\sigma = " << sigma << '$' << endl;
//  int nevt = lround( data->sumEntries() );
//  fsig = frac2.getVal();
//  float esig = frac2.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << width << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
  if ( psig ) model.plotOn( mframe, Components( "crybf" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) model.plotOn( mframe, Components( "polyn" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  model.plotOn( mframe );
  mframe->Draw();
  cout << mframe->GetYaxis()->GetXmin() << " <---> " << mframe->GetYaxis()->GetXmax() << endl;
  cout << mframe->GetMinimum() << " <---> " << mframe->GetMaximum() << endl;
/*
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
//  h->Draw( "same" );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
  if ( wwid ) printSigma( mframe, sigma, xmin, xmax, xfac, yfac );
  printInfo( mframe, lv, xmin, xmax, xfac, yfac );
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

}

void fitCB2Pol2( const string& head,
                 const string& name,
                 const string& pset, int year,
                 const string& dset, float lumi,
                 const string& ptit,
                 const string& atit,
                 const vector<const InfoLabel*>& lv,
                 int nbin, double xmin, double xmax, float xfac, float yfac,
                 double mass, double sig1, double sig2,
                 double tcut, double tpow, double fsig,
                 bool wwid = true, bool pbkg = true, bool psig = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar rmass( "rmass", "Resonance  Mass" , mass, xmin, xmax, "GeV" );
  RooRealVar widt1( "widt1", "Resonance Widt1" , sig1, sig1/10, sig1*2, "GeV" );
  RooRealVar widt2( "widt2", "Resonance Widt2" , sig2, sig2/10, sig2*2, "GeV" );
  RooRealVar sfcut( "sfcut", "Cut Distance"    , tcut, 0.0, tcut*2, "GeV" );
  RooRealVar power( "power", "Tail Power"      , tpow, 0.0, tpow*10, "." );
  RooRealVar coef1( "coef1", "Coefficicent 1"  , 0.0, -100.0, 100.0, "GeV^-1" );
  RooRealVar coef2( "coef2", "Coefficicent 1"  , 0.0,  -10.0,  10.0, "GeV^-2" );
  RooRealVar frac1( "frac1", "signal fraction1",  0.5, 0.0, 1.0 );
//  RooRealVar frac2( "frac2", "signal fraction2", fsig, 0.0, 1.0 );
  // create the fit model components: CB and polynomial PDFs
  RooCBShape cryb1( "cryb1", "Crystal Ball Lineshape1", m, rmass, widt1,
                                                           sfcut, power );
  RooCBShape cryb2( "cryb2", "Crystal Ball Lineshape2", m, rmass, widt2,
                                                           sfcut, power );
  RooPolynomial polyn( "polyn", "Polynomial Distribution", m, RooArgList( coef1, coef2 ) );
  // combine them using addition
  RooAddPdf cbsum( "cbsum", "Signal", RooArgList( cryb1, cryb2 ), frac1 );
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, m );
  int nevr = data->numEntries();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs / 10.0, nevs * 10.0 );
  RooExtendPdf exsig( "exsig", "esig", cbsum, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", polyn, nbkg );
  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbkg ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( cbsum, polyn ), frac2 );
  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data, Extended( kTRUE ) );
//  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
/*
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "widt1: " << widt1.getVal() << " +/- " << widt1.getError() << endl;
  cout << "widt2: " << widt2.getVal() << " +/- " << widt2.getError() << endl;
  cout << "sfcut: " << sfcut.getVal() << " +/- " << sfcut.getError() << endl;
  cout << "power: " << power.getVal() << " +/- " << power.getError() << endl;
  cout << "coef1: " << coef1.getVal() << " +/- " << coef1.getError() << endl;
  cout << "coef2: " << coef2.getVal() << " +/- " << coef2.getError() << endl;
  cout << "frac1: " << frac1.getVal() << " +/- " << frac1.getError() << endl;
//  cout << "frac2: " << frac2.getVal() << " +/- " << frac2.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
*/
  float f1 = frac1.getVal();
  float sigma = sqrt( ( pow( widt1.getVal(), 2 ) *         f1   ) +
                      ( pow( widt2.getVal(), 2 ) * ( 1.0 - f1 ) ) );
  ofstream os( ( nfit + ".tex" ).c_str() );
  print( rmass, os );
  print( widt1, os );
  print( widt2, os );
  print( sfcut, os );
  print( power, os );
  print( coef1, os );
  print( coef2, os );
  print( frac1, os );
  print( nsig , os );
  print( nbkg , os );
  print( nsig , os, 1.0 / lumi );
  cout << " sigma = " << sigma << endl;
  os << "$\\sigma = " << sigma << '$' << endl;
//  int nevt = lround( data->sumEntries() );
//  fsig = frac2.getVal();
//  float esig = frac2.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << sigma << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
  if ( psig ) model.plotOn( mframe, Components( "cbsum" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) model.plotOn( mframe, Components( "polyn" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  model.plotOn( mframe );
  mframe->Draw();
  cout << mframe->GetYaxis()->GetXmin() << " <---> " << mframe->GetYaxis()->GetXmax() << endl;
  cout << mframe->GetMinimum() << " <---> " << mframe->GetMaximum() << endl;
/*
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
//  h->Draw( "same" );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
  if ( wwid ) printSigma( mframe, sigma, xmin, xmax, xfac, yfac );
  printInfo( mframe, lv, xmin, xmax, xfac, yfac );
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

}

void fitCBChP2 ( const string& head,
                 const string& name,
                 const string& pset, int year,
                 const string& dset, float lumi,
                 const string& ptit,
                 const string& atit,
                 const vector<const InfoLabel*>& lv,
                 int nbin, double xmin, double xmax, float xfac, float yfac,
                 double mass, double sigm, double tcut, double tpow,
                 double fsig,
                 bool wwid = true, bool pbkg = true, bool psig = false,
                 bool line = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar rmass( "rmass", "Resonance  Mass" , mass, xmin, xmax, "GeV" );
  RooRealVar width( "width", "Resonance Width" , sigm, sigm/10, sigm*5, "GeV" );
  RooRealVar sfcut( "sfcut", "Cut Distance"    , tcut, 0.0, tcut*10, "GeV" );
  RooRealVar power( "power", "Tail Power"      , tpow, 0.0, tpow*10, "." );
  RooRealVar coef1( "coef1", "Coefficicent 1"  , 0.0, -100.0, 100.0, "GeV^-1" );
  RooRealVar coef2( "coef2", "Coefficicent 1"  , 0.0,  -10.0,  10.0, "GeV^-2" );
//  RooRealVar frac1( "frac1", "signal fraction1",  0.5, 0.0, 1.0 );
//  RooRealVar frac2( "frac2", "signal fraction2", fsig, 0.0, 1.0 );
  // create the fit model components: CB and polynomial PDFs
  RooCBShape crybf( "crybf", "Crystal Ball Lineshape", m, rmass, width,
                                                          sfcut, power );
  RooChebychev polyn( "polyn", "Polynomial Distribution", m, RooArgList( coef1, coef2 ) );
  // combine them using addition
//  RooAddPdf cbsum( "cbsum", "Signal", RooArgList( cryb1, cryb2 ), frac1 );
//  RooAddPdf tmpdf( "tmpdf", "Signal + Background tmp", RooArgList( cbsum, polyn ), frac2 );
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, m );
  int nevr = data->numEntries();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
/*
  tmpdf.fitTo( *data );
  fsig = frac2.getVal();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
  cout << nevs << ' ' << nevb << endl;
//  RooAddPdf& model = tmpdf;
  float errs = 10 * sqrt( nevs );
  float errb = 10 * sqrt( nevb );
  if ( errs < errb ) errs = errb;
  else               errb = errs;
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs - errs, nevs + errs );
  RooExtendPdf exsig( "exsig", "esig", cbsum, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb - errb, nevb + errb );
  RooExtendPdf exbkg( "exbkg", "ebkg", polyn, nbkg );
*/
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs / 10.0, nevs * 10.0 );
  RooExtendPdf exsig( "exsig", "esig", crybf, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", polyn, nbkg );
  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbkg ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( cbsum, polyn ), frac2 );
  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data, Extended( kTRUE ) );
//  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
/*
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "width: " << width.getVal() << " +/- " << width.getError() << endl;
  cout << "sfcut: " << sfcut.getVal() << " +/- " << sfcut.getError() << endl;
  cout << "power: " << power.getVal() << " +/- " << power.getError() << endl;
  cout << "coef1: " << coef1.getVal() << " +/- " << coef1.getError() << endl;
  cout << "coef2: " << coef2.getVal() << " +/- " << coef2.getError() << endl;
//  cout << "frac1: " << frac1.getVal() << " +/- " << frac1.getError() << endl;
//  cout << "frac2: " << frac2.getVal() << " +/- " << frac2.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
*/
  float sigma = width.getVal();
  ofstream os( ( nfit + ".tex" ).c_str() );
  print( rmass, os );
  print( width, os );
  print( sfcut, os );
  print( power, os );
  print( coef1, os );
  print( coef2, os );
//  print( frac1, os );
  print( nsig , os );
  print( nbkg , os );
  print( nsig , os, 1.0 / lumi );
//  float f1 = frac1.getVal();
//  float sigma = sqrt( ( pow( widt1.getVal(), 2 ) *         f1   ) +
//                      ( pow( widt2.getVal(), 2 ) * ( 1.0 - f1 ) ) );
  os << "$\\sigma = " << sigma << '$' << endl;
  cout << " sigma = " << sigma << endl;
//  int nevt = lround( data->sumEntries() );
//  fsig = frac2.getVal();
//  float esig = frac2.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << width << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
  if ( psig ) model.plotOn( mframe, Components( "crybf" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) model.plotOn( mframe, Components( "polyn" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  model.plotOn( mframe );
  mframe->Draw();
  cout << mframe->GetYaxis()->GetXmin() << " <---> " << mframe->GetYaxis()->GetXmax() << endl;
  cout << mframe->GetMinimum() << " <---> " << mframe->GetMaximum() << endl;
/*
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
//  h->Draw( "same" );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
  if ( wwid ) printSigma( mframe, sigma, xmin, xmax, xfac, yfac );
  printInfo( mframe, lv, xmin, xmax, xfac, yfac );
  yfac -= ( 0.07 + ( 0.06 * lv.size() ) );
  if ( line ) {
    printLine( mframe, xmin, xmax, xfac, yfac, "Total Fit" );
    if ( psig ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                           "Signal"           , kDashed, csig );
    if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                           "Combinatorial bkg", kDashed, cbkg );
  }
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

}

void fitCBGChP2 ( const string& head,
                  const string& name,
                  const string& pset, int year,
                  const string& dset, float lumi,
                  const string& ptit,
                  const string& atit,
                  const vector<const InfoLabel*>& lv,
                  const string& cbfl, const string& gaul,
                  int nbin, double xmin, double xmax, float xfac, float yfac,
                  double mcbf, double scbf, double tcut, double tpow,
                  double mgau, double sgau,
                  double fcbf, double fgau,
//                  bool wwid = true,
                  bool pbkg = true, bool pscb = false, bool pgau = false,
                  bool line = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar cmass( "cmass", "CBFun  Mass" , mcbf, xmin, xmax, "GeV" );
  RooRealVar width( "width", "CBFun Width" , scbf, scbf/10, scbf*5, "GeV" );
  RooRealVar gmass( "gmass", "Gauss  Mass" , mgau, xmin, xmax, "GeV" );
  RooRealVar sigma( "sigma", "Gauss Width" , sgau, sgau/10, sgau*5, "GeV" );
  RooRealVar sfcut( "sfcut", "Cut Distance"    , tcut, 0.0, tcut*10, "GeV" );
  RooRealVar power( "power", "Tail Power"      , tpow, 0.0, tpow*10, "." );
  RooRealVar coef1( "coef1", "Coefficicent 1"  , 0.0, -100.0, 100.0, "GeV^-1" );
  RooRealVar coef2( "coef2", "Coefficicent 1"  , 0.0,  -10.0,  10.0, "GeV^-2" );
//  RooRealVar frac1( "frac1", "signal fraction1",  0.5, 0.0, 1.0 );
//  RooRealVar frac2( "frac2", "signal fraction2", fsig, 0.0, 1.0 );
  // create the fit model components: CB and polynomial PDFs
  RooCBShape crybf( "crybf", "Crystal Ball Lineshape", m, cmass, width,
                                                          sfcut, power );
  RooGaussian gauss( "gauss", "Gaussian Distribution", m, gmass, sigma );
  RooChebychev polyn( "polyn", "Polynomial Distribution", m, RooArgList( coef1, coef2 ) );
  // combine them using addition
//  RooAddPdf cbsum( "cbsum", "Signal", RooArgList( cryb1, cryb2 ), frac1 );
//  RooAddPdf tmpdf( "tmpdf", "Signal + Background tmp", RooArgList( cbsum, polyn ), frac2 );
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, m );
  int nevr = data->numEntries();
  float nevc = nevr * fcbf;
  float nevg = nevr * fgau;
  float nevb = nevr * ( 1.0 - ( fcbf + fgau ) );
  RooRealVar ncbf( "ncbf", "nscrybf", nevc, nevc / 10.0, nevc * 10.0 );
  RooExtendPdf excbf( "excbf", "ecbf", crybf, ncbf );
  RooRealVar ngau( "ngau", "nsgauss", nevg, nevg / 10.0, nevg * 10.0 );
  RooExtendPdf exgau( "exgau", "egau", gauss, ngau );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", polyn, nbkg );
  RooAddPdf model( "model", "Signal + Background", RooArgList( excbf, exgau, exbkg ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( cbsum, polyn ), frac2 );
  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data, Extended( kTRUE ) );
//  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
/*
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "width: " << width.getVal() << " +/- " << width.getError() << endl;
  cout << "sfcut: " << sfcut.getVal() << " +/- " << sfcut.getError() << endl;
  cout << "power: " << power.getVal() << " +/- " << power.getError() << endl;
  cout << "coef1: " << coef1.getVal() << " +/- " << coef1.getError() << endl;
  cout << "coef2: " << coef2.getVal() << " +/- " << coef2.getError() << endl;
//  cout << "frac1: " << frac1.getVal() << " +/- " << frac1.getError() << endl;
//  cout << "frac2: " << frac2.getVal() << " +/- " << frac2.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
*/
//  float sigma = width.getVal();
  ofstream os( ( nfit + ".tex" ).c_str() );
  print( cmass, os );
  print( width, os );
  print( sfcut, os );
  print( power, os );
  print( gmass, os );
  print( sigma, os );
  print( coef1, os );
  print( coef2, os );
//  print( frac1, os );
  print( ncbf , os );
  print( ngau , os );
  print( nbkg , os );
  print( ncbf , os, 1.0 / lumi );
  print( ngau , os, 1.0 / lumi );
//  float f1 = frac1.getVal();
//  float sigma = sqrt( ( pow( widt1.getVal(), 2 ) *         f1   ) +
//                      ( pow( widt2.getVal(), 2 ) * ( 1.0 - f1 ) ) );
//  os << "$\\sigma = " << sigma << '$' << endl;
//  cout << " sigma = " << sigma << endl;
//  int nevt = lround( data->sumEntries() );
//  fsig = frac2.getVal();
//  float esig = frac2.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << width << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
  if ( pscb ) model.plotOn( mframe, Components( "crybf" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pgau ) model.plotOn( mframe, Components( "gauss" ),
                            LineStyle( kDashed ), LineColor( csi2 ) );
  if ( pbkg ) model.plotOn( mframe, Components( "polyn" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  model.plotOn( mframe );
  mframe->Draw();
  cout << mframe->GetYaxis()->GetXmin() << " <---> " << mframe->GetYaxis()->GetXmax() << endl;
  cout << mframe->GetMinimum() << " <---> " << mframe->GetMaximum() << endl;
/*
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
//  h->Draw( "same" );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
//  if ( wwid ) printSigma( mframe, sigma, xmin, xmax, xfac, yfac );
  printInfo( mframe, lv, xmin, xmax, xfac, yfac );
  yfac -= ( 0.07 + ( 0.06 * lv.size() ) );
  if ( line ) {
    printLine( mframe, xmin, xmax, xfac, yfac, "Total Fit" );
    if ( pscb ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                           cbfl               , kDashed, csig );
    if ( pgau ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                           gaul               , kDashed, csi2 );
    if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                           "Combinatorial bkg", kDashed, cbkg );
  }
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

}

void fitCB2ChP2( const string& head,
                 const string& name,
                 const string& pset, int year,
                 const string& dset, float lumi,
                 const string& ptit,
                 const string& atit,
                 const vector<const InfoLabel*>& lv,
                 int nbin, double xmin, double xmax, float xfac, float yfac,
                 double mass, double sig1, double sig2,
                 double tcut, double tpow, double fsig,
                 bool wwid = true, bool pbkg = true, bool psig = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar rmass( "rmass", "Resonance  Mass" , mass, xmin, xmax, "GeV" );
  RooRealVar widt1( "widt1", "Resonance Widt1" , sig1, sig1/10, sig1*5, "GeV" );
  RooRealVar widt2( "widt2", "Resonance Widt2" , sig2, sig2/10, sig2*5, "GeV" );
  RooRealVar sfcut( "sfcut", "Cut Distance"    , tcut, 0.0, tcut*10, "GeV" );
  RooRealVar power( "power", "Tail Power"      , tpow, 0.0, tpow*10, "." );
  RooRealVar coef1( "coef1", "Coefficicent 1"  , 0.0, -100.0, 100.0, "GeV^-1" );
  RooRealVar coef2( "coef2", "Coefficicent 1"  , 0.0,  -10.0,  10.0, "GeV^-2" );
  RooRealVar frac1( "frac1", "signal fraction1",  0.5, 0.0, 1.0 );
//  RooRealVar frac2( "frac2", "signal fraction2", fsig, 0.0, 1.0 );
  // create the fit model components: CB and polynomial PDFs
  RooCBShape cryb1( "cryb1", "Crystal Ball Lineshape1", m, rmass, widt1,
                                                           sfcut, power );
  RooCBShape cryb2( "cryb2", "Crystal Ball Lineshape2", m, rmass, widt2,
                                                           sfcut, power );
  RooChebychev polyn( "polyn", "Polynomial Distribution", m, RooArgList( coef1, coef2 ) );
  // combine them using addition
  RooAddPdf cbsum( "cbsum", "Signal", RooArgList( cryb1, cryb2 ), frac1 );
//  RooAddPdf tmpdf( "tmpdf", "Signal + Background tmp", RooArgList( cbsum, polyn ), frac2 );
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, m );
  int nevr = data->numEntries();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
/*
  tmpdf.fitTo( *data );
  fsig = frac2.getVal();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
  cout << nevs << ' ' << nevb << endl;
//  RooAddPdf& model = tmpdf;
  float errs = 10 * sqrt( nevs );
  float errb = 10 * sqrt( nevb );
  if ( errs < errb ) errs = errb;
  else               errb = errs;
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs - errs, nevs + errs );
  RooExtendPdf exsig( "exsig", "esig", cbsum, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb - errb, nevb + errb );
  RooExtendPdf exbkg( "exbkg", "ebkg", polyn, nbkg );
*/
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs / 10.0, nevs * 10.0 );
  RooExtendPdf exsig( "exsig", "esig", cbsum, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", polyn, nbkg );
  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbkg ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( cbsum, polyn ), frac2 );
  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data, Extended( kTRUE ) );
//  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
/*
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "widt1: " << widt1.getVal() << " +/- " << widt1.getError() << endl;
  cout << "widt2: " << widt2.getVal() << " +/- " << widt2.getError() << endl;
  cout << "sfcut: " << sfcut.getVal() << " +/- " << sfcut.getError() << endl;
  cout << "power: " << power.getVal() << " +/- " << power.getError() << endl;
  cout << "coef1: " << coef1.getVal() << " +/- " << coef1.getError() << endl;
  cout << "coef2: " << coef2.getVal() << " +/- " << coef2.getError() << endl;
  cout << "frac1: " << frac1.getVal() << " +/- " << frac1.getError() << endl;
//  cout << "frac2: " << frac2.getVal() << " +/- " << frac2.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
*/
  ofstream os( ( nfit + ".tex" ).c_str() );
  print( rmass, os );
  print( widt1, os );
  print( widt2, os );
  print( sfcut, os );
  print( power, os );
  print( coef1, os );
  print( coef2, os );
  print( frac1, os );
  print( nsig , os );
  print( nbkg , os );
  print( nsig , os, 1.0 / lumi );
  float f1 = frac1.getVal();
  float sigma = sqrt( ( pow( widt1.getVal(), 2 ) *         f1   ) +
                      ( pow( widt2.getVal(), 2 ) * ( 1.0 - f1 ) ) );
  os << "$\\sigma = " << sigma << '$' << endl;
  cout << " sigma = " << sigma << endl;
//  int nevt = lround( data->sumEntries() );
//  fsig = frac2.getVal();
//  float esig = frac2.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << sigma << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
  if ( psig ) model.plotOn( mframe, Components( "cbsum" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) model.plotOn( mframe, Components( "polyn" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  model.plotOn( mframe );
  mframe->Draw();
  cout << mframe->GetYaxis()->GetXmin() << " <---> " << mframe->GetYaxis()->GetXmax() << endl;
  cout << mframe->GetMinimum() << " <---> " << mframe->GetMaximum() << endl;
/*
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
//  h->Draw( "same" );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
  if ( wwid ) printSigma( mframe, sigma, xmin, xmax, xfac, yfac );
  printInfo( mframe, lv, xmin, xmax, xfac, yfac );
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

}

void fitCB4ChP2( const string& head,
                 const string& name,
                 const string& pset, int year,
                 const string& dset, float lumi,
                 const string& ptit,
                 const string& atit,
                 const vector<const InfoLabel*>& lv,
                 int nbin, double xmin, double xmax, float xfac, float yfac,
                 double mas1, double siga, double sigb,
                 double cut1, double pow1,
                 double mas2, double sig2, double cut2, double pow2,
                 double mas3, double sig3, double cut3, double pow3,
                 double fsig,
                 bool wwid = true, bool pbkg = true, bool psig = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar rmas1( "rmas1", "Resonance  Mass1", mas1, xmin, xmax, "GeV" );
  RooRealVar rmas2( "rmas2", "Resonance  Mass2", mas2, xmin, xmax, "GeV" );
  RooRealVar rmas3( "rmas3", "Resonance  Mass3", mas3, xmin, xmax, "GeV" );
  RooRealVar widta( "widta", "Resonance WidtA" , siga, siga/10, siga*2, "GeV" );
  RooRealVar widtb( "widtb", "Resonance WidtB" , sigb, sigb/10, sigb*2, "GeV" );
  RooRealVar widt2( "widt2", "Resonance Widt2" , sig2, sig2/10, sig2*2, "GeV" );
  RooRealVar widt3( "widt3", "Resonance Widt3" , sig3, sig3/10, sig3*2, "GeV" );
  RooRealVar scut1( "scut1", "Cut Distance1"   , cut1, 0.0, cut1*5, "GeV" );
  RooRealVar powt1( "powt1", "Tail Power1"     , pow1, 0.0, pow1*10, "." );
  RooRealVar scut2( "scut2", "Cut Distance2"   , cut2, 0.0, cut2*5, "GeV" );
  RooRealVar powt2( "powt2", "Tail Power2"     , pow2, 0.0, pow2*10, "." );
  RooRealVar scut3( "scut3", "Cut Distance3"   , cut3, 0.0, cut3*5, "GeV" );
  RooRealVar powt3( "powt3", "Tail Power3"     , pow3, 0.0, pow3*10, "." );
  RooRealVar coef1( "coef1", "Coefficicent 1"  , 0.0, -100.0, 100.0, "GeV^-1" );
  RooRealVar coef2( "coef2", "Coefficicent 1"  , 0.0,  -10.0,  10.0, "GeV^-2" );
  RooRealVar fraca( "fraca", "signal fraction1",  0.5, 0.0, 1.0 );
//  RooRealVar frac1( "frac1", "signal fraction1",  0.5, 0.0, 1.0 );
//  RooRealVar frac2( "frac2", "signal fraction2",  0.3, 0.0, 1.0 );
//  RooRealVar fracs( "fracs", "signal fractions", fsig, 0.0, 1.0 );
  // create the fit model components: CB and polynomial PDFs
  RooCBShape cryba( "cryba", "Crystal Ball Lineshapea", m, rmas1, widta,
                                                           scut1, powt1 );
  RooCBShape crybb( "crybb", "Crystal Ball Lineshapeb", m, rmas1, widtb,
                                                           scut1, powt1 );
  RooCBShape cryb2( "cryb2", "Crystal Ball Lineshape2", m, rmas2, widt2,
                                                           scut2, powt2 );
  RooCBShape cryb3( "cryb3", "Crystal Ball Lineshape3", m, rmas3, widt3,
                                                           scut3, powt3 );
  RooChebychev polyn( "polyn", "Polynomial Distribution", m, RooArgList( coef1, coef2 ) );
  // combine them using addition
  RooAddPdf cbsum( "cbsum", "Signal1", RooArgList( cryba, crybb ), fraca );
//  RooAddPdf cbtot( "cbtot", "SignalT", RooArgList( cbsum, cryb2, cryb3 ), RooArgList( frac1, frac2 ) );
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, m );
  int nevr = data->numEntries();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
  float nev1 = nevs * 0.5;
  float nev2 = nevs * 0.3;
  float nev3 = nevs * 0.2;
//  RooRealVar nsig( "nsig", "nsignal", nevs, nevs / 10.0, nevs * 10.0 );
//  RooExtendPdf exsig( "exsig", "esig", cbtot, nsig );
  RooRealVar nsi1( "nsi1", "nsignal1", nev1, nev1 / 10.0, nev1 * 10.0 );
  RooExtendPdf exsi1( "exsi1", "esig1", cbsum, nsi1 );
  RooRealVar nsi2( "nsi2", "nsignal2", nev2, nev2 / 10.0, nev2 * 10.0 );
  RooExtendPdf exsi2( "exsi2", "esig2", cryb2, nsi2 );
  RooRealVar nsi3( "nsi3", "nsignal3", nev3, nev3 / 10.0, nev3 * 10.0 );
  RooExtendPdf exsi3( "exsi3", "esig3", cryb3, nsi3 );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", polyn, nbkg );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbkg ) );
  RooAddPdf model( "model", "Signal + Background", RooArgList( exsi1, exsi2, exsi3, exbkg ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( cbtot, polyn ), fracs );
  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data, Extended( kTRUE ) );
//  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
/*
  cout << "rmas1: " << rmas1.getVal() << " +/- " << rmas1.getError() << endl;
  cout << "rmas2: " << rmas2.getVal() << " +/- " << rmas2.getError() << endl;
  cout << "rmas3: " << rmas3.getVal() << " +/- " << rmas3.getError() << endl;
  cout << "widta: " << widta.getVal() << " +/- " << widta.getError() << endl;
  cout << "widtb: " << widtb.getVal() << " +/- " << widtb.getError() << endl;
  cout << "widt2: " << widt2.getVal() << " +/- " << widt2.getError() << endl;
  cout << "widt3: " << widt3.getVal() << " +/- " << widt3.getError() << endl;
  cout << "scut1: " << scut1.getVal() << " +/- " << scut1.getError() << endl;
  cout << "powt1: " << powt1.getVal() << " +/- " << powt1.getError() << endl;
  cout << "scut2: " << scut2.getVal() << " +/- " << scut2.getError() << endl;
  cout << "powt2: " << powt2.getVal() << " +/- " << powt2.getError() << endl;
  cout << "scut3: " << scut3.getVal() << " +/- " << scut3.getError() << endl;
  cout << "powt3: " << powt3.getVal() << " +/- " << powt3.getError() << endl;
  cout << "coef1: " << coef1.getVal() << " +/- " << coef1.getError() << endl;
  cout << "coef2: " << coef2.getVal() << " +/- " << coef2.getError() << endl;
  cout << "fraca: " << fraca.getVal() << " +/- " << fraca.getError() << endl;
//  cout << "frac1: " << frac1.getVal() << " +/- " << frac1.getError() << endl;
//  cout << "frac2: " << frac2.getVal() << " +/- " << frac2.getError() << endl;
//  cout << "fracs: " << fracs.getVal() << " +/- " << fracs.getError() << endl;
//  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nsi1:  " << nsi1 .getVal() << " +/- " << nsi1 .getError() << endl;
  cout << "nsi2:  " << nsi2 .getVal() << " +/- " << nsi2 .getError() << endl;
  cout << "nsi3:  " << nsi3 .getVal() << " +/- " << nsi3 .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
//  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
//                    << nsig .getError() / lumi << endl;
  cout << "nrm1:  " << nsi1 .getVal() / lumi << " +/- "
                    << nsi1 .getError() / lumi << endl;
  cout << "nrm2:  " << nsi2 .getVal() / lumi << " +/- "
                    << nsi2 .getError() / lumi << endl;
  cout << "nrm3:  " << nsi3 .getVal() / lumi << " +/- "
                    << nsi3 .getError() / lumi << endl;
*/
  float fa = fraca.getVal();
  float sigma = sqrt( ( pow( widta.getVal(), 2 ) *         fa   ) +
                      ( pow( widtb.getVal(), 2 ) * ( 1.0 - fa ) ) );
  ofstream os( ( nfit + ".tex" ).c_str() );
  print( rmas1, os );
  print( rmas2, os );
  print( rmas3, os );
  print( widta, os );
  print( widtb, os );
  print( widt2, os );
  print( widt3, os );
  print( scut1, os );
  print( powt1, os );
  print( scut2, os );
  print( powt2, os );
  print( scut3, os );
  print( powt3, os );
  print( coef1, os );
  print( coef2, os );
  print( fraca, os );
  print( nsi1 , os );
  print( nsi2 , os );
  print( nsi3 , os );
  print( nbkg , os );
  print( nsi1 , os, 1.0 / lumi );
  print( nsi2 , os, 1.0 / lumi );
  print( nsi3 , os, 1.0 / lumi );
  cout << " sigma = " << sigma << endl;
  os << "$\\sigma = " << sigma << '$' << endl;
//  int nevt = lround( data->sumEntries() );
//  fsig = frac2.getVal();
//  float esig = frac2.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << sigma << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
  if ( psig ) model.plotOn( mframe, Components( "cbsum" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) model.plotOn( mframe, Components( "polyn" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  model.plotOn( mframe );
  mframe->Draw();
/*
  cout << mframe->GetYaxis()->GetXmin() << " <---> " << mframe->GetYaxis()->GetXmax() << endl;
  cout << mframe->GetMinimum() << " <---> " << mframe->GetMaximum() << endl;
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
//  h->Draw( "same" );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
  if ( wwid ) printSigma( mframe, sigma, xmin, xmax, xfac, yfac );
  printInfo( mframe, lv, xmin, xmax, xfac, yfac );
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

}

void fitGauss2Exp( const string& head,
                   const string& name,
                   const string& pset, int year,
                   const string& dset, float lumi,
                   const string& ptit,
                   const string& atit,
                   const vector<const InfoLabel*>& lv,
                   int nbin, double xmin, double xmax, float xfac, float yfac,
                   double mass, double sig1, double sig2,
                   double expc, double fsig,
                   bool pbkg = true, bool psig = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar rmass( "rmass", "Resonance  Mass", mass, xmin, xmax, "GeV" );
  RooRealVar widt1( "widt1", "Resonance Widt1", sig1, sig1/3, sig1*3, "GeV" );
  RooRealVar widt2( "widt2", "Resonance Widt2", sig2, sig2/3, sig2*3, "GeV" );
  RooRealVar bgpar( "bgpar", "Background Par" , expc, -100.0, 0.0, "GeV" );
  RooRealVar frac1( "frac1", "signal fraction1",  0.5, 0.0, 1.0 );
//  RooRealVar frac2( "frac2", "signal fraction2", fsig, 0.0, 1.0 );
  // create the fit model components: Gaussian and exponential PDFs
  RooGaussian gauss1( "gauss1", "Signal Distribution1", m, rmass, widt1 );
  RooGaussian gauss2( "gauss2", "Signal Distribution2", m, rmass, widt2 );
  RooExponential expbg( "expbg", "Background Distribution", m, bgpar );
  // combine them using addition
  RooAddPdf gsum2( "gsum2", "Signal", RooArgList( gauss1, gauss2 ), frac1 );
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, m );
  int nevr = data->numEntries();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs / 10.0, nevs * 10.0 );
  RooExtendPdf exsig( "exsig", "esig", gsum2, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", expbg, nbkg );
  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbkg ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( gsum2, expbg ), frac2 );
//  RooAddPdf model( "model", "Signal + Background",
//                   RooArgList( gauss1, gauss2, expbg ),
//                   RooArgList( frac1, frac2 ) );
  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data, Extended( kTRUE ) );
//  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
/*
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "widt1: " << widt1.getVal() << " +/- " << widt1.getError() << endl;
  cout << "widt2: " << widt2.getVal() << " +/- " << widt2.getError() << endl;
  cout << "bgpar: " << bgpar.getVal() << " +/- " << bgpar.getError() << endl;
  cout << "frac1: " << frac1.getVal() << " +/- " << frac1.getError() << endl;
//  cout << "frac2: " << frac2.getVal() << " +/- " << frac2.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
*/
  float f1 = frac1.getVal();
  float sigma = sqrt( ( pow( widt1.getVal(), 2 ) *         f1   ) +
                      ( pow( widt2.getVal(), 2 ) * ( 1.0 - f1 ) ) );
  ofstream os( ( nfit + ".tex" ).c_str() );
  print( rmass, os );
  print( widt1, os );
  print( widt2, os );
  print( bgpar, os );
  print( frac1, os );
  print( nsig , os );
  print( nbkg , os );
  print( nsig , os, 1.0 / lumi );
  cout << " sigma = " << sigma << endl;
  os << "$\\sigma = " << sigma << '$' << endl;
//  int nevt = lround( data->sumEntries() );
//  fsig = frac2.getVal();
//  float esig = frac2.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << sigma << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
  if ( psig ) model.plotOn( mframe, Components( "gsum2" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) model.plotOn( mframe, Components( "expbg" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  model.plotOn( mframe );
  mframe->Draw();
/*
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
  printInfo( mframe, lv, xmin, xmax, xfac, yfac );
//  float xfac = 0.65;
//  float yfac = 0.78 - ( 0.06 * lv.size() );
  yfac -= ( 0.07 + ( 0.06 * lv.size() ) );
  printLine( mframe, xmin, xmax, xfac, yfac, "Total Fit" );
  if ( psig ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Signal"           , kDashed, csig );
  if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Combinatorial bkg", kDashed, cbkg );
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

}

void fitGauss2ExpE( const string& head,
                    const string& name,
                    const string& pset, int year,
                    const string& dset, float lumi,
                    const string& ptit,
                    const string& atit,
                    const vector<const InfoLabel*>& lv, const string& misrl,
                    int nbin, double xmin, double xmax, float xfac, float yfac,
                    double mass, double sig1, double sig2,
                    double expc, double mbkg, double sbkg, double fbkg,
                    double fsig,
                    bool pbkg = true, bool psig = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar rmass( "rmass", "Resonance  Mass", mass, xmin, xmax, "GeV" );
  RooRealVar widt1( "widt1", "Resonance Widt1", sig1, 0.003, sig1*3, "GeV" );
  RooRealVar widt2( "widt2", "Resonance Widt2", sig2, 0.003, sig2*3, "GeV" );
  RooRealVar bgpar( "bgpar", "Background Par" , expc, -100.0, 0.0, "GeV" );
  RooRealVar xmass( "xmass", "BG XMass", mbkg, mbkg-0.1, mbkg+0.1, "GeV" );
  RooRealVar xwidt( "xwidt", "BG XWidt", sbkg, sbkg/5, sbkg*5, "GeV" );
  RooRealVar frac1( "frac1", "signal fraction1",  0.5, 0.0, 1.0 );
  RooRealVar fracb( "fracb", "backgr fraction" , fbkg, 0.0, 1.0 );
//  RooRealVar frac2( "frac2", "signal fraction2", fsig, 0.0, 1.0 );
  // create the fit model components: Gaussian and exponential PDFs
  RooGaussian gauss1( "gauss1", "Signal Distribution1", m, rmass, widt1 );
  RooGaussian gauss2( "gauss2", "Signal Distribution2", m, rmass, widt2 );
  RooExponential expbg( "expbg", "Background Combinatorial", m, bgpar );
//  RooGaussian gausbg( "gausbg", "Background Partial Reco", m, xmass, xwidt );
//  RooAbsPdf* efptr = bindPdf( "erfun", ROOT::Math::Erf, m, xmass, xwidt );
//  TF1* tferf = new TF1( "tferf", "1-erf([1]*(x-[0]))", xmin, xmax );
//  RooAbsPdf* efptr = bindPdf( "erfun", tferf, m, xmass, xwidt );
//  ErrF errFunct( xmin, xmax );
//  RooAbsPdf* efptr = bindPdf( "erfun", [errFunct]( double x, double threshold, double slope ) { return errFunct( x, threshold, slope );}, m, xmass, xwidt );
  RooAbsPdf* efptr = bindPdf( "erfun", errFunct, m, xmass, xwidt );
  RooAbsPdf& erfun = *efptr;
  // combine them using addition
  RooAddPdf gsum2( "gsum2", "Signal", RooArgList( gauss1, gauss2 ), frac1 );
  RooAddPdf bgsum( "bgsum", "Background", RooArgList( expbg, erfun ), fracb );
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, m );
  int nevr = data->numEntries();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs / 10.0, nevs * 10.0 );
  RooExtendPdf exsig( "exsig", "esig", gsum2, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", bgsum, nbkg );
  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbkg ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( gsum2, bgsum ), frac2 );
  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data, Extended( kTRUE ) );
//  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
/*
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "widt1: " << widt1.getVal() << " +/- " << widt1.getError() << endl;
  cout << "widt2: " << widt2.getVal() << " +/- " << widt2.getError() << endl;
  cout << "bgpar: " << bgpar.getVal() << " +/- " << bgpar.getError() << endl;
  cout << "xmass: " << xmass.getVal() << " +/- " << xmass.getError() << endl;
  cout << "xwidt: " << xwidt.getVal() << " +/- " << xwidt.getError() << endl;
  cout << "frac1: " << frac1.getVal() << " +/- " << frac1.getError() << endl;
  cout << "fracb: " << fracb.getVal() << " +/- " << fracb.getError() << endl;
//  cout << "frac2: " << frac2.getVal() << " +/- " << frac2.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
*/
  float f1 = frac1.getVal();
  float sigma = sqrt( ( pow( widt1.getVal(), 2 ) *         f1   ) +
                      ( pow( widt2.getVal(), 2 ) * ( 1.0 - f1 ) ) );
  ofstream os( ( nfit + ".tex" ).c_str() );
  print( rmass, os );
  print( widt1, os );
  print( widt2, os );
  print( bgpar, os );
  print( xmass, os );
  print( xwidt, os );
  print( frac1, os );
  print( fracb, os );
  print( nsig , os );
  print( nbkg , os );
  print( nsig , os, 1.0 / lumi );
  cout << " sigma = " << sigma << endl;
  os << "$\\sigma = " << sigma << '$' << endl;
//  int nevt = lround( data->sumEntries() );
//  fsig = frac2.getVal();
//  float esig = frac2.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << sigma << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
  if ( psig ) model.plotOn( mframe, Components( "gsum2" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) model.plotOn( mframe, Components( "expbg" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) model.plotOn( mframe, Components( "erfun" ),
                            LineStyle( kDashed ), LineColor( cmis ) );
  model.plotOn( mframe );
  mframe->Draw();
/*
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
  printInfo( mframe, lv, xmin, xmax, xfac, yfac );
//  float xfac = 0.65;
//  float yfac = 0.78 - ( 0.06 * lv.size() );
  yfac -= ( 0.07 + ( 0.06 * lv.size() ) );
  printLine( mframe, xmin, xmax, xfac, yfac, "Total Fit" );
  if ( psig ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Signal"           , kDashed, csig );
  if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Combinatorial bkg", kDashed, cbkg );
  if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         misrl, kDashed, cmis );
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

}

void fitGauss2ExpE2D( const string& head,
                      const string& name,
                      const string& pset, int year,
                      const string& dset, float lumi,
                      const string& ptit,
                      const string& atit,
                      const vector<const InfoLabel*>& lv, const string& misrl,
                      int nbin, double xmin, double xmax,
                      float xfac, float yfac, float xfct, float yfct,
                      double mass, double sig1, double sig2,
                      double expc, double mbkg, double sbkg, double fbkg,
                      double fsig,
                      bool pbkg = true, bool psig = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  float dmin = -0.05;
  float dmax =  0.25;
  int dbin = 60;
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar d( "d", "Flight Distance", dmin, dmax, "cm" );
  RooRealVar e( "e", "Distance Error", 0.00001, 0.1, "cm" );
//  RooRealVar e( "e", "Distance Error", 0.0, 0.25, "cm" );
  RooRealVar rmass( "rmass", "Resonance  Mass", mass, xmin, xmax, "GeV" );
  RooRealVar widt1( "widt1", "Resonance Widt1", sig1, 0.001, sig1*10, "GeV" );
  RooRealVar widt2( "widt2", "Resonance Widt2", sig2, 0.001, sig2*10, "GeV" );
  RooRealVar bgpar( "bgpar", "Background Par" , expc, -100.0, 0.0, "GeV" );
  RooRealVar xmass( "xmass", "BG XMass", mbkg, mbkg-0.1, mbkg+0.1, "GeV" );
  RooRealVar xwidt( "xwidt", "BG XWidt", sbkg, sbkg/5, sbkg*5, "GeV" );
  RooRealVar frac1( "frac1", "signal fraction1",  0.5, 0.0, 1.0 );
  RooRealVar fracb( "fracb", "backgr fraction" , fbkg, 0.0, 1.0 );
//  RooRealVar frac2( "frac2", "signal fraction2", fsig, 0.0, 1.0 );
  // create the fit model components: Gaussian and exponential PDFs
  RooGaussian gauss1( "gauss1", "Signal Distribution1", m, rmass, widt1 );
  RooGaussian gauss2( "gauss2", "Signal Distribution2", m, rmass, widt2 );
  RooExponential expbg( "expbg", "Background Combinatorial", m, bgpar );
//  RooGaussian gausbg( "gausbg", "Background Partial Reco", m, xmass, xwidt );
//  RooAbsPdf* efptr = bindPdf( "erfun", ROOT::Math::Erf, m, xmass, xwidt );
//  TF1* tferf = new TF1( "tferf", "1-erf([1]*(x-[0]))", xmin, xmax );
//  RooAbsPdf* efptr = bindPdf( "erfun", tferf, m, xmass, xwidt );
//  ErrF errFunct( xmin, xmax );
//  RooAbsPdf* efptr = bindPdf( "erfun", [errFunct]( double x, double threshold, double slope ) { return errFunct( x, threshold, slope );}, m, xmass, xwidt );
  RooAbsPdf* efptr = bindPdf( "erfun", errFunct, m, xmass, xwidt );
  RooAbsPdf& erfun = *efptr;
  // combine them using addition
  RooAddPdf gsum2( "gsum2", "Signal", RooArgList( gauss1, gauss2 ), frac1 );
  RooAddPdf bgsum( "bgsum", "Background", RooArgList( expbg, erfun ), fracb );
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, RooArgList( m, d, e ) );
  int nevr = data->numEntries();
  float nevs = nevr * fsig;
  float nevb = nevr - nevs;
  float nevc = nevb * fbkg;
  float nevf = nevb - nevc;
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs / 10.0, nevs * 10.0 );
  RooExtendPdf exsig( "exsig", "esig", gsum2, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", bgsum, nbkg );
  RooRealVar nbgf( "nbgf", "nbackgf", nevf, nevf / 10.0, nevf * 10.0 );
  RooExtendPdf exbgf( "exbgf", "ebgf", erfun, nbgf );
  RooRealVar nbgc( "nbgc", "nbackgc", nevc, nevc / 10.0, nevc * 10.0 );
  RooExtendPdf exbgc( "exbgc", "ebgc", expbg, nbgc );
  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbgf, exbgc ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbkg ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( gsum2, bgsum ), frac2 );
  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data, Extended( kTRUE ) );
//  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "widt1: " << widt1.getVal() << " +/- " << widt1.getError() << endl;
  cout << "widt2: " << widt2.getVal() << " +/- " << widt2.getError() << endl;
  cout << "bgpar: " << bgpar.getVal() << " +/- " << bgpar.getError() << endl;
  cout << "xmass: " << xmass.getVal() << " +/- " << xmass.getError() << endl;
  cout << "xwidt: " << xwidt.getVal() << " +/- " << xwidt.getError() << endl;
  cout << "frac1: " << frac1.getVal() << " +/- " << frac1.getError() << endl;
  cout << "fracb: " << fracb.getVal() << " +/- " << fracb.getError() << endl;
//  cout << "frac2: " << frac2.getVal() << " +/- " << frac2.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "nbgf:  " << nbgf .getVal() << " +/- " << nbgf .getError() << endl;
  cout << "nbgc:  " << nbgc .getVal() << " +/- " << nbgc .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
  float f1 = frac1.getVal();
  float sigma = sqrt( ( pow( widt1.getVal(), 2 ) *         f1   ) +
                      ( pow( widt2.getVal(), 2 ) * ( 1.0 - f1 ) ) );
  cout << " sigma = " << sigma << endl;


  TCanvas* tmp = new TCanvas( "tmp", nfit.c_str(), 800, 600 );
  tmp->SetLeftMargin( 0.15 );
  tmp->cd();
  m.setBins( nbin );
  RooPlot* mfrtmp = m.frame();
  mfrtmp->SetTitle( ptit.c_str() );
  mfrtmp->SetXTitle( atit.c_str() );
  data->plotOn( mfrtmp );
  if ( psig ) model.plotOn( mfrtmp, Components( "gsum2" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) model.plotOn( mfrtmp, Components( "expbg" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) model.plotOn( mfrtmp, Components( "erfun" ),
                            LineStyle( kDashed ), LineColor( cmis ) );
  model.plotOn( mfrtmp );
  mfrtmp->Draw();
//  return;


//  int i;
//  cin >> i;
//  cout << i << endl;
//  RooRealVar btcom( "btcom", "Combin Bg Life", -100.0, -1000.0, 0.0, "cm" );
//  RooRealVar btfly( "btfly", "Flying Bg Life", -20.0, -500.0, 0.0, "cm" );
//  RooRealVar stime( "stime", "Signal Life"   , -20.0, -500.0, 0.0, "cm" );
  RooRealVar btcom( "btcom", "Combin Bg Life", 0.01, 0.0, 0.1, "cm" );
  RooRealVar bqcom( "bqcom", "Combin BR Life", 0.02, 0.0, 0.1, "cm" );
  RooRealVar brcom( "brcom", "Combin BR Life", 0.03, 0.0, 0.1, "cm" );
  RooRealVar blcom( "blcom", "Combin BL Life", 0.02, 0.0, 0.1, "cm" );
//  btcom.setConstant( kTRUE );
  RooRealVar btfly( "btfly", "Flying Bg Life", 0.05, 0.0, 0.15, "cm" );
  RooRealVar stime( "stime", "Signal Life"   , 0.05, 0.0, 0.15, "cm" );
  RooRealVar dmean( "dmean", "Distance Res Mean", 0.0, 0.0, 0.01, "cm" );
  dmean.setConstant( kTRUE );
  RooRealVar dreso( "dreso", "Distance Resolution", 0.002, 0.00001, 0.1, "cm" );
  RooRealVar drfac( "drfac", "Distance Resolution", 1.0, 0.05, 20.0, "cm" );
  RooRealVar widtp( "widtp", "Prompt width", 0.005, 0.00001, 0.1, "cm" );
//  widtp.setConstant( kTRUE );
  RooRealVar fracp( "fracp", "Prompt Core Fraction", 0.5, 0.0, 1.0, "." );
  RooRealVar fracq( "fracq", "Q-Tail Fraction", 0.5, 0.0, 1.0, "." );
  RooRealVar fracr( "fracr", "R-Tail Fraction", 0.7, 0.0, 1.0, "." );
  RooRealVar fract( "fract", "PTails Fraction", 0.1, 0.0, 1.0, "." );
/*
  RooGaussian gaussr( "gaussr", "Distance Resolution", d, dmean, dreso );
  RooExponential etbgc( "etbgc", "Combin Background Decay", d, btcom );
  RooExponential etbgf( "etbgf", "Flying Background Decay", d, btfly );
  RooExponential etsig( "etsig", "Signal Decay", d, stime );
  RooNumConvPdf ecbgc( "ecbgc", "Combin Convolved Decay", d, etbgc, gaussr );
  RooNumConvPdf ecbgf( "ecbgf", "Flying Convolved Decay", d, etbgf, gaussr );
  RooNumConvPdf ecsig( "ecsig", "Signal Convolved Decay", d, etsig, gaussr );
*/
//  RooGaussModel gaussr( "gaussr", "Distance Resolution", d, dmean, drfac, e );
  RooGaussModel gaussr( "gaussr", "Distance Resolution", d, dmean, dreso );
  RooDecay ecsig( "ecsig", "Signal Convolved Decay", d, stime, gaussr, RooDecay::SingleSided );
  RooDecay ecbgf( "ecbgf", "Flying Convolved Decay", d, btfly, gaussr, RooDecay::SingleSided );
  RooDecay ecbgr( "ecbgr", "Combin Convolved Decay", d, brcom, gaussr, RooDecay::SingleSided );
  RooDecay ecbgq( "ecbgq", "Combin Convolved Decay", d, bqcom, gaussr, RooDecay::SingleSided );
  RooDecay ecbgl( "ecbgl", "Combin Convolved Decay", d, blcom, gaussr, RooDecay::Flipped );
  RooDecay ecbgd( "ecbgd", "Combin Convolved Decay", d, btcom, gaussr, RooDecay::DoubleSided );
  RooGaussian gaussc( "gaussc", "Prompt Core Distribution", d, dmean, dreso );
  RooGaussian gausst( "gausst", "Prompt Tail Distribution", d, dmean, widtp );
  RooAddPdf gaussp( "gaussp", "Convolved Tails sum", RooArgList( gaussc, gausst ), fracp );

  RooAddPdf ecbgs( "ecbgs", "Convolved Tails sum", RooArgList( ecbgr, ecbgq ), fracq );
//  RooAddPdf ecbgt( "ecbgt", "Convolved Tails sum", RooArgList( ecbgs, ecbgl ), fracr );
  RooAddPdf ecbgt( "ecbgt", "Convolved Tails sum", RooArgList( ecbgr, ecbgl ), fracr );
  RooAddPdf ddbgr( "ddbgr", "Combin Total Decay", RooArgList( ecbgr, gaussp ), fract );
  RooAddPdf ddbgl( "ddbgl", "Combin Total Decay", RooArgList( ecbgl, gaussp ), fract );
  RooAddPdf ddbgt( "ddbgt", "Combin Total Decay", RooArgList( ecbgt, gaussp ), fract );
  RooAddPdf ddbgd( "ddbgd", "Combin Total Decay", RooArgList( ecbgd, gaussp ), fract );

  RooProdPdf mtsig( "mtsig", "Signal 2D", RooArgSet( gsum2, ecsig ) );
  RooProdPdf mtbgf( "mtbgf", "Flying 2D", RooArgSet( erfun, ecbgf ) );
  RooProdPdf mtbgr( "mtbgr", "Combin 2D", RooArgSet( expbg, ddbgr ) );
  RooProdPdf mtbgl( "mtbgl", "Combin 2D", RooArgSet( expbg, ddbgl ) );
  RooProdPdf mtbgt( "mtbgt", "Combin 2D", RooArgSet( expbg, ddbgt ) );
  RooProdPdf mtbgd( "mtbgd", "Combin 2D", RooArgSet( expbg, ddbgd ) );
  RooAddPdf mtbkg( "mtbkg", "Backgr 2D", RooArgList( mtbgd, mtbgf ), fracb );

  RooExtendPdf edsig( "edsig", "edsig", mtsig, nsig );
  RooExtendPdf edbgf( "edbgf", "edbgf", mtbgf, nbgf );
  RooExtendPdf edbkg( "edbkg", "edbkg", mtbkg, nbkg );
  RooExtendPdf edbgr( "edbgr", "edbgr", mtbgr, nbgc );
  RooExtendPdf edbgl( "edbgl", "edbgl", mtbgl, nbgc );
  RooExtendPdf edbgt( "edbgt", "edbgt", mtbgt, nbgc );
  RooExtendPdf edbgd( "edbgd", "edbgd", mtbgd, nbgc );
  RooAddPdf mdisr( "mdisr", "Signal + Background", RooArgList( edsig, edbgf, edbgr ) );
  RooAddPdf mdisl( "mdisl", "Signal + Background", RooArgList( edsig, edbgf, edbgl ) );
  RooAddPdf mdiss( "mdiss", "Signal + Background", RooArgList( edsig, edbgf, edbgt ) );
  RooAddPdf mdisd( "mdisd", "Signal + Background", RooArgList( edsig, edbgf, edbgd ) );
//  RooAddPdf mdist( "mdist", "Signal + Background", RooArgList( edsig, edbkg ) );

  TString sdcut = "m>5.5";
//  TString sdcut = "m>5.22 && m<5.32";
//  TString sdcut = "m>";
//  sdcut += rmass.getVal() + ( 5 * sigma );
  RooAbsData* side = data->reduce( sdcut );

//  RooAbsData* dtmp = data;
//  RooAddPdf& mdist = mdisd;
//  RooAddPdf& mdtmp = mdisl;
//  Bool_t extm = kTRUE;
  RooAbsData* dtmp = side;
  RooAddPdf& mdist = mdiss;
  RooAddPdf& mdtmp = ddbgt;
  Bool_t extm = kFALSE;

  rmass.setConstant( kTRUE );
  widt1.setConstant( kTRUE );
  widt2.setConstant( kTRUE );
  frac1.setConstant( kTRUE );
  bgpar.setConstant( kTRUE );
  xmass.setConstant( kTRUE );
  xwidt.setConstant( kTRUE );
//  fracb.setConstant( kTRUE );
//  nsig .setConstant( kTRUE );
//  nbkg .setConstant( kTRUE );
//  nbgf .setConstant( kTRUE );
//  nbgc .setConstant( kTRUE );

//  fracp.setConstant( kTRUE );
//  fracr.setConstant( kTRUE );
//  fract.setConstant( kTRUE );
//  brcom.setConstant( kTRUE );
//  blcom.setConstant( kTRUE );
//  btcom.setConstant( kTRUE );

//  ddbgt.fitTo( *side );

//  mdist.fitTo( *data, ConditionalObservables( e ), Extended( kTRUE ) );
//  mdtmp.fitTo( *data, Extended( kTRUE ) );
  mdtmp.fitTo( *dtmp, extm );

  TCanvas* t2d = new TCanvas( "t2d", nfit.c_str(), 800, 600 );
  t2d->SetLeftMargin( 0.15 );
  t2d->cd();
  t2d->SetLogy();
  d.setBins( dbin );
  RooPlot* dfrtm2 = d.frame();
  dfrtm2->SetXTitle( atit.c_str() );
//  m.setBins( nbin );
//  RooPlot* dfrtm2 = m.frame();
//  dfrtm2->SetXTitle( "ct [cm]" );
  dfrtm2->SetTitle( ptit.c_str() );
  dtmp->plotOn( dfrtm2 );
/*
  if ( psig ) mdtmp.plotOn( dfrtm2, Components( "mtsig" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) mdtmp.plotOn( dfrtm2, Components( "mtbgl" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) mdtmp.plotOn( dfrtm2, Components( "mtbgr" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) mdtmp.plotOn( dfrtm2, Components( "mtbgt" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) mdtmp.plotOn( dfrtm2, Components( "mtbgd" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) mdtmp.plotOn( dfrtm2, Components( "mtbgf" ),
                            LineStyle( kDashed ), LineColor( cmis ) );
*/
  // single expo
  mdtmp.plotOn( dfrtm2, Components( "gaussp" ),
                        LineStyle( kDashed ), LineColor( csig ) );
  mdtmp.plotOn( dfrtm2, Components( "ecbgr" ),
                        LineStyle( kDashed ), LineColor( cbkg ) );
  mdtmp.plotOn( dfrtm2, Components( "ecbgl" ),
                        LineStyle( kDashed ), LineColor( cmis ) );
  mdtmp.plotOn( dfrtm2 );
  dfrtm2->Draw();
  dfrtm2->SetMaximum( 1.2 * dfrtm2->GetMaximum() );
  dfrtm2->SetMinimum( 0.5 );

  ofstream ot( ( nfit + "tmpFit.tex" ).c_str() );
  ot << "% " << sdcut << endl;
  ot << "% " << data->numEntries() << endl;
  ot << "% " << dtmp->numEntries() << endl;
  print( rmass, ot );
  print( widt1, ot );
  print( widt2, ot );
  print( bgpar, ot );
  print( xmass, ot );
  print( xwidt, ot );
  print( frac1, ot );
  print( fracb, ot );
  print( brcom, ot );
  print( bqcom, ot );
  print( blcom, ot );
  print( btcom, ot );
  print( btfly, ot );
  print( stime, ot );
  print( dmean, ot );
  print( dreso, ot );
  print( widtp, ot );
  print( drfac, ot );
//  print( widtp, ot );
  print( fracp, ot );
  print( fracr, ot );
  print( fracq, ot );
  print( fract, ot );
  print( nsig , ot );
  print( nbkg , ot );
  print( nbgf , ot );
  print( nbgc , ot );
  print( nsig , ot, 1.0 / lumi );
//  nsig .setConstant( kTRUE );
//  nbgf .setConstant( kTRUE );
//  nbgc .setConstant( kTRUE );
//  double feg = fract.getVal();
//  fract.setVal( 2 * feg / ( 1 + feg ) );
//  btcom.setVal( blcom.getVal() );
  fracp.setConstant( kTRUE );
  fracr.setConstant( kTRUE );
  fracq.setConstant( kTRUE );
  fract.setConstant( kTRUE );
  brcom.setConstant( kTRUE );
  bqcom.setConstant( kTRUE );
  blcom.setConstant( kTRUE );
  btcom.setConstant( kTRUE );
  widtp.setConstant( kTRUE );
  dreso.setConstant( kTRUE );
//  btfly.setConstant( kTRUE );
//  stime.setConstant( kTRUE );

//  return;
  mdist.fitTo( *data, Extended( kTRUE ) );


  string nf2d = nfit + "Dist";
  cout << nfit << " results" << endl;
/*
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "widt1: " << widt1.getVal() << " +/- " << widt1.getError() << endl;
  cout << "widt2: " << widt2.getVal() << " +/- " << widt2.getError() << endl;
  cout << "bgpar: " << bgpar.getVal() << " +/- " << bgpar.getError() << endl;
  cout << "xmass: " << xmass.getVal() << " +/- " << xmass.getError() << endl;
  cout << "xwidt: " << xwidt.getVal() << " +/- " << xwidt.getError() << endl;
  cout << "frac1: " << frac1.getVal() << " +/- " << frac1.getError() << endl;
  cout << "fracb: " << fracb.getVal() << " +/- " << fracb.getError() << endl;
//  cout << "frac2: " << frac2.getVal() << " +/- " << frac2.getError() << endl;
  cout << "btcom: " << btcom.getVal() << " +/- " << btcom.getError() << endl;
  cout << "btfly: " << btfly.getVal() << " +/- " << btfly.getError() << endl;
  cout << "stime: " << stime.getVal() << " +/- " << stime.getError() << endl;
  cout << "dmean: " << dmean.getVal() << " +/- " << dmean.getError() << endl;
  cout << "dreso: " << dreso.getVal() << " +/- " << dreso.getError() << endl;
  cout << "widtp: " << widtp.getVal() << " +/- " << widtp.getError() << endl;
  cout << "fract: " << fract.getVal() << " +/- " << fract.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
*/
  f1 = frac1.getVal();
  sigma = sqrt( ( pow( widt1.getVal(), 2 ) *         f1   ) +
                ( pow( widt2.getVal(), 2 ) * ( 1.0 - f1 ) ) );
  ofstream os( ( nfit + ".tex" ).c_str() );
  cout << sdcut << endl;
  os << "% " << sdcut << endl;
  os << "% " << data->numEntries() << endl;
  os << "% " << dtmp->numEntries() << endl;
  print( rmass, os );
  print( widt1, os );
  print( widt2, os );
  print( bgpar, os );
  print( xmass, os );
  print( xwidt, os );
  print( frac1, os );
  print( fracb, os );
  print( brcom, os );
  print( bqcom, os );
  print( blcom, os );
  print( btcom, os );
  print( btfly, os );
  print( stime, os );
  print( dmean, os );
  print( dreso, os );
  print( widtp, os );
  print( drfac, os );
//  print( widtp, os );
  print( fracp, os );
  print( fracr, os );
  print( fracq, os );
  print( fract, os );
  print( nsig , os );
  print( nbkg , os );
  print( nbgf , os );
  print( nbgc , os );
  print( nsig , os, 1.0 / lumi );
  cout << " sigma = " << sigma << endl;
  os << "$\\sigma = " << sigma << '$' << endl;

//  int nevt = lround( data->sumEntries() );
//  fsig = frac2.getVal();
//  float esig = frac2.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << sigma << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
//  if ( psig ) model.plotOn( mframe, Components( "gsum2" ),
//                            LineStyle( kDashed ), LineColor( csig ) );
//  if ( pbkg ) model.plotOn( mframe, Components( "expbg" ),
//                            LineStyle( kDashed ), LineColor( cbkg ) );
//  if ( pbkg ) model.plotOn( mframe, Components( "erfun" ),
//                            LineStyle( kDashed ), LineColor( cmis ) );
//  model.plotOn( mframe );
  if ( psig ) mdist.plotOn( mframe, Components( "mtsig" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) mdist.plotOn( mframe, Components( "mtbgr" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) mdist.plotOn( mframe, Components( "mtbgl" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) mdist.plotOn( mframe, Components( "mtbgt" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) mdist.plotOn( mframe, Components( "mtbgd" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) mdist.plotOn( mframe, Components( "mtbgf" ),
                            LineStyle( kDashed ), LineColor( cmis ) );
  mdist.plotOn( mframe );
  mframe->Draw();
/*
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
  printInfo( mframe, lv, xmin, xmax, xfac, yfac );
//  float xfac = 0.65;
//  float yfac = 0.78 - ( 0.06 * lv.size() );
  yfac -= ( 0.07 + ( 0.06 * lv.size() ) );
  printLine( mframe, xmin, xmax, xfac, yfac, "Total Fit" );
  if ( psig ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Signal"           , kDashed, csig );
  if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Combinatorial bkg", kDashed, cbkg );
  if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         misrl, kDashed, cmis );
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

  TCanvas* c2d = new TCanvas( nf2d.c_str(), nf2d.c_str(), 800, 600 );
  c2d->SetLeftMargin( 0.15 );
  c2d->cd();
  c2d->SetLogy();
  d.setBins( dbin );
  RooPlot* dframe = d.frame();
  dframe->SetTitle( ptit.c_str() );
  dframe->SetXTitle( "ct [cm]" );
  data->plotOn( dframe );
  if ( psig ) mdist.plotOn( dframe, Components( "mtsig" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) mdist.plotOn( dframe, Components( "mtbgr" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) mdist.plotOn( dframe, Components( "mtbgl" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) mdist.plotOn( dframe, Components( "mtbgt" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) mdist.plotOn( dframe, Components( "mtbgd" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) mdist.plotOn( dframe, Components( "mtbgf" ),
                            LineStyle( kDashed ), LineColor( cmis ) );
  mdist.plotOn( dframe );
  dframe->Draw();
/*
  float ymax = 1.2 * dframe->GetMaximum();
  cout << "DFRAME: " << ymax << endl;
  dframe->SetMaximum( ymax );
  dframe->SetMinimum( ymax / 2500.0 );
//  t.DrawLatex( ( dmin * 0.45 ) + ( dmax * 0.55 ), dframe->GetMaximum() * 1.01, s.c_str() );
  float ef = 1.01;
//  t.DrawLatex( ( dmin * 0.45 ) + ( dmax * 0.55 ), 10.0, s.c_str() );
//  t.DrawLatex( ( dmin * 0.45 ) + ( dmax * 0.55 ), pow( dframe->GetMaximum(), ef ), s.c_str() );
  t.DrawLatex( ( dmin * 0.45 ) + ( dmax * 0.55 ), pow( dframe->GetMinimum(), ef - 1.0 ) * pow( dframe->GetMaximum(), ef ), s.c_str() );
*/
  float topMargin = 1.2;
  float range = 2500.0;
  float ymax = topMargin * dframe->GetMaximum();
  dframe->SetMaximum( ymax );
  dframe->SetMinimum( ymax / range );
  printTitleLog( dframe, year, lumi, dmin, dmax );
  xfac = xfct;
  yfac = yfct;
  printInfo( dframe, lv, dmin, dmax, xfac, yfac, InfoLabel::log );
  yfac -= ( 0.07 + ( 0.06 * lv.size() ) );
  printLine( dframe, dmin, dmax, xfac, yfac, "Total Fit",
             kSolid, kBlue, InfoLabel::log );
  if ( psig ) printLine( dframe, dmin, dmax, xfac, yfac -= 0.04,
                         "Signal",
                         kDashed, csig, InfoLabel::log );
  if ( pbkg ) printLine( dframe, dmin, dmax, xfac, yfac -= 0.04,
                         "Combinatorial bkg",
                         kDashed, cbkg, InfoLabel::log );
  if ( pbkg ) printLine( dframe, dmin, dmax, xfac, yfac -= 0.04,
                         misrl,
                         kDashed, cmis, InfoLabel::log );
//  nf2d += ".pdf";
  c2d->Print( ( nf2d + ".pdf" ).c_str() );
  c2d->Print( ( nf2d + ".png" ).c_str() );

}

void fitGauss2ExpG( const string& head,
                    const string& name,
                    const string& pset, int year,
                    const string& dset, float lumi,
                    const string& ptit,
                    const string& atit,
                    const vector<const InfoLabel*>& lv, const string& misrl,
                    int nbin, double xmin, double xmax, float xfac, float yfac,
                    double mass, double sig1, double sig2,
                    double expc, double mbkg, double sbkg, double fbkg,
                    double fsig,
                    bool pbkg = true, bool psig = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar rmass( "rmass", "Resonance  Mass", mass, xmin, xmax, "GeV" );
  RooRealVar widt1( "widt1", "Resonance Widt1", sig1, 0.003, sig1*3, "GeV" );
  RooRealVar widt2( "widt2", "Resonance Widt2", sig2, 0.003, sig2*3, "GeV" );
  RooRealVar bgpar( "bgpar", "Background Par" , expc, -100.0, 0.0, "GeV" );
  RooRealVar xmass( "xmass", "BG XMass", mbkg, mbkg-0.1, mbkg+0.1, "GeV" );
  RooRealVar xwidt( "xwidt", "BG XWidt", sbkg, sbkg/5, sbkg*5, "GeV" );
  RooRealVar frac1( "frac1", "signal fraction1",  0.5, 0.0, 1.0 );
  RooRealVar fracb( "fracb", "backgr fraction" , fbkg, 0.0, 1.0 );
//  RooRealVar frac2( "frac2", "signal fraction2", fsig, 0.0, 1.0 );
  // create the fit model components: Gaussian and exponential PDFs
  RooGaussian gauss1( "gauss1", "Signal Distribution1", m, rmass, widt1 );
  RooGaussian gauss2( "gauss2", "Signal Distribution2", m, rmass, widt2 );
  RooExponential expbg( "expbg", "Background Combinatorial", m, bgpar );
  RooGaussian gausbg( "gausbg", "Background Partial Reco", m, xmass, xwidt );
  // combine them using addition
  RooAddPdf gsum2( "gsum2", "Signal", RooArgList( gauss1, gauss2 ), frac1 );
  RooAddPdf bgsum( "bgsum", "Background", RooArgList( expbg, gausbg ), fracb );
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, m );
  int nevr = data->numEntries();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs / 10.0, nevs * 10.0 );
  RooExtendPdf exsig( "exsig", "esig", gsum2, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", bgsum, nbkg );
  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbkg ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( gsum2, bgsum ), frac2 );
  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data, Extended( kTRUE ) );
//  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
/*
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "widt1: " << widt1.getVal() << " +/- " << widt1.getError() << endl;
  cout << "widt2: " << widt2.getVal() << " +/- " << widt2.getError() << endl;
  cout << "bgpar: " << bgpar.getVal() << " +/- " << bgpar.getError() << endl;
  cout << "xmass: " << xmass.getVal() << " +/- " << xmass.getError() << endl;
  cout << "xwidt: " << xwidt.getVal() << " +/- " << xwidt.getError() << endl;
  cout << "frac1: " << frac1.getVal() << " +/- " << frac1.getError() << endl;
  cout << "fracb: " << fracb.getVal() << " +/- " << fracb.getError() << endl;
//  cout << "frac2: " << frac2.getVal() << " +/- " << frac2.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
*/
  float f1 = frac1.getVal();
  float sigma = sqrt( ( pow( widt1.getVal(), 2 ) *         f1   ) +
                      ( pow( widt2.getVal(), 2 ) * ( 1.0 - f1 ) ) );
  ofstream os( ( nfit + ".tex" ).c_str() );
  print( rmass, os );
  print( widt1, os );
  print( widt2, os );
  print( bgpar, os );
  print( xmass, os );
  print( xwidt, os );
  print( frac1, os );
  print( fracb, os );
  print( nsig , os );
  print( nbkg , os );
  print( nsig , os, 1.0 / lumi );
  cout << " sigma = " << sigma << endl;
  os << "$\\sigma = " << sigma << '$' << endl;
//  int nevt = lround( data->sumEntries() );
//  fsig = frac2.getVal();
//  float esig = frac2.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << sigma << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
  if ( psig ) model.plotOn( mframe, Components( "gsum2" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) model.plotOn( mframe, Components( "expbg" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) model.plotOn( mframe, Components( "gausbg" ),
                            LineStyle( kDashed ), LineColor( cmis ) );
  model.plotOn( mframe );
  mframe->Draw();
/*
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
  printInfo( mframe, lv, xmin, xmax, xfac, yfac );
//  float xfac = 0.65;
//  float yfac = 0.78 - ( 0.06 * lv.size() );
  yfac -= ( 0.07 + ( 0.06 * lv.size() ) );
  printLine( mframe, xmin, xmax, xfac, yfac, "Total Fit" );
  if ( psig ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Signal"           , kDashed, csig );
  if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Combinatorial bkg", kDashed, cbkg );
  if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         misrl, kDashed, cmis );
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

}

void fitGaussExp( const string& head,
                  const string& name,
                  const string& pset, int year,
                  const string& dset, float lumi,
                  const string& ptit,
                  const string& atit,
                  const vector<const InfoLabel*>& lv,
                  int nbin, double xmin, double xmax, float xfac, float yfac,
                  double mass, double ssig,
                  double expc, double fsig,
                  bool pbkg = true, bool psig = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
//  cout << file << endl;
//  return;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar rmass( "rmass", "Resonance  Mass", mass, xmin, xmax, "GeV" );
  RooRealVar width( "width", "Resonance Widts", ssig, ssig/3, ssig*3, "GeV" );
  RooRealVar bgpar( "bgpar", "Background Par" , expc, -100.0, 0.0, "GeV" );
  RooRealVar fracs( "fracs", "signal fraction", fsig, 0.0, 1.0 );
  // create the fit model components: Gaussian and exponential PDFs
  RooGaussian gauss( "gauss", "Signal Distribution1", m, rmass, width );
  RooExponential expbg( "expbg", "Background Distribution", m, bgpar );
  // combine them using addition
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, m );
  int nevr = data->numEntries();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs / 10.0, nevs * 10.0 );
  RooExtendPdf exsig( "exsig", "esig", gauss, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", expbg, nbkg );
  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbkg ) );
/*
*/
//  RooAddPdf model( "model", "Signal + Background", RooArgList( gauss, expbg ), fracs );
//  RooAddPdf model( "model", "Signal + Background",
//                   RooArgList( gauss1, gauss2, expbg ),
//                   RooArgList( frac1, frac2 ) );

  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
/*
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "width: " << width.getVal() << " +/- " << width.getError() << endl;
  cout << "bgpar: " << bgpar.getVal() << " +/- " << bgpar.getError() << endl;
//  cout << "fracs: " << fracs.getVal() << " +/- " << fracs.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
*/
  float sigma = width.getVal();
  ofstream os( ( nfit + ".tex" ).c_str() );
  print( rmass, os );
  print( width, os );
  print( bgpar, os );
  print( fracs , os );
  print( fracs , os, nevr );
  print( fracs , os, nevr / lumi);
  print( nsig , os );
  print( nbkg , os );
  print( nsig , os, 1.0 / lumi );
  cout << " sigma = " << sigma << endl;
  os << "$\\sigma = " << sigma << '$' << endl;
  os << "$\\sigma = " << sigma << '$' << endl;
//  int nevt = lround( data->sumEntries() );
//  fsig = fracs.getVal();
//  float esig = fracs.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << sigma << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
  if ( psig ) model.plotOn( mframe, Components( "gauss" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) model.plotOn( mframe, Components( "expbg" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  model.plotOn( mframe, LineColor( kBlue ) );
  mframe->Draw();
/*
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
  printInfo( mframe, lv, xmin, xmax, xfac, yfac );
//  float 
  yfac -= ( 0.07 + ( 0.06 * lv.size() ) );
  printLine( mframe, xmin, xmax, xfac, yfac, "Total Fit" );
  if ( psig ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Signal"           , kDashed, csig );
  if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Combinatorial bkg", kDashed, cbkg );
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

}

void fitGaussExpE( const string& head,
                   const string& name,
                   const string& pset, int year,
                   const string& dset, float lumi,
                   const string& ptit,
                   const string& atit,
                   const vector<const InfoLabel*>& lv, const string& misrl,
                   int nbin, double xmin, double xmax, float xfac, float yfac,
                   double mass, double ssig,
                   double expc, double mbkg, double sbkg, double fbkg,
                   double fsig,
                   bool pbkg = true, bool psig = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar rmass( "rmass", "Resonance  Mass", mass, xmin, xmax, "GeV" );
  RooRealVar width( "width", "Resonance Width", ssig, ssig/3, ssig*3, "GeV" );
  RooRealVar bgpar( "bgpar", "Background Par" , expc, -100.0, 0.0, "GeV" );
  RooRealVar xmass( "xmass", "BG XMass", mbkg, mbkg-0.1, mbkg+0.1, "GeV" );
  RooRealVar xwidt( "xwidt", "BG XWidt", sbkg, sbkg/5, sbkg*5, "GeV" );
  RooRealVar fracb( "fracb", "backgr fraction", fbkg, 0.0, 1.0 );
//  RooRealVar fracs( "fracs", "signal fraction", fsig, 0.0, 1.0 );
  // create the fit model components: Gaussian and exponential PDFs
  RooGaussian gausig( "gausig", "Signal Distribution1", m, rmass, width );
  RooExponential expbg( "expbg", "Background Combinatorial", m, bgpar );
//  RooGaussian gausbg( "gausbg", "Background Partial Reco", m, xmass, xwidt );
  RooAbsPdf* efptr = bindPdf( "erfun", errFunct, m, xmass, xwidt );
  RooAbsPdf& erfun = *efptr;
  // combine them using addition
  RooAddPdf bgsum( "bgsum", "Background", RooArgList( expbg, erfun ), fracb );
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, m );
  int nevr = data->numEntries();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs / 10.0, nevs * 10.0 );
  RooExtendPdf exsig( "exsig", "esig", gausig, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", bgsum, nbkg );
  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbkg ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( gausig, bgsum ), fracs );
  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data, Extended( kTRUE ) );
//  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
/*
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "width: " << width.getVal() << " +/- " << width.getError() << endl;
  cout << "bgpar: " << bgpar.getVal() << " +/- " << bgpar.getError() << endl;
  cout << "xmass: " << xmass.getVal() << " +/- " << xmass.getError() << endl;
  cout << "xwidt: " << xwidt.getVal() << " +/- " << xwidt.getError() << endl;
  cout << "fracb: " << fracb.getVal() << " +/- " << fracb.getError() << endl;
//  cout << "fracs: " << fracs.getVal() << " +/- " << fracs.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
*/
  ofstream os( ( nfit + ".tex" ).c_str() );
  print( rmass, os );
  print( width, os );
  print( bgpar, os );
  print( xmass, os );
  print( xwidt, os );
  print( fracb, os );
  print( nsig , os );
  print( nbkg , os );
  print( nsig , os, 1.0 / lumi );
  float sigma = width.getVal();
  cout << " sigma = " << sigma << endl;
  os << "$\\sigma = " << sigma << '$' << endl;
//  int nevt = lround( data->sumEntries() );
//  fsig = fracs.getVal();
//  float esig = fracs.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << sigma << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
  if ( psig ) model.plotOn( mframe, Components( "gsum2" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) model.plotOn( mframe, Components( "expbg" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) model.plotOn( mframe, Components( "erfun" ),
                            LineStyle( kDashed ), LineColor( cmis ) );
  model.plotOn( mframe );
  mframe->Draw();
/*
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
  printInfo( mframe, lv, xmin, xmax, xfac, yfac );
//  float xfac = 0.65;
//  float yfac = 0.78 - ( 0.06 * lv.size() );
  yfac -= ( 0.07 + ( 0.06 * lv.size() ) );
  printLine( mframe, xmin, xmax, xfac, yfac, "Total Fit" );
  if ( psig ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Signal"           , kDashed, csig );
  if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Combinatorial bkg", kDashed, cbkg );
  if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         misrl, kDashed, cmis );
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

}

void fitGaussExpG( const string& head,
                   const string& name,
                   const string& pset, int year,
                   const string& dset, float lumi,
                   const string& ptit,
                   const string& atit,
                   const vector<const InfoLabel*>& lv, const string& misrl,
                   int nbin, double xmin, double xmax, float xfac, float yfac,
                   double mass, double ssig,
                   double expc, double mbkg, double sbkg, double fbkg,
                   double fsig,
                   bool pbkg = true, bool psig = false ) {
  string nfit = name + pset + to_string( year ) + dset;
  string file = head + nfit;
  const char* list = file.c_str();
  // define the data variables and fit model parameters
  RooRealVar m( "m", "Reconstructed Mass", xmin, xmax, "GeV" );
  RooRealVar rmass( "rmass", "Resonance  Mass", mass, xmin, xmax, "GeV" );
  RooRealVar width( "width", "Resonance Width", ssig, ssig/3, ssig*3, "GeV" );
  RooRealVar bgpar( "bgpar", "Background Par" , expc, -100.0, 0.0, "GeV" );
  RooRealVar xmass( "xmass", "BG XMass", mbkg, mbkg-0.1, mbkg+0.1, "GeV" );
  RooRealVar xwidt( "xwidt", "BG XWidt", sbkg, sbkg/5, sbkg*5, "GeV" );
  RooRealVar fracb( "fracb", "backgr fraction", fbkg, 0.0, 1.0 );
//  RooRealVar fracs( "fracs", "signal fraction", fsig, 0.0, 1.0 );
  // create the fit model components: Gaussian and exponential PDFs
  RooGaussian gausig( "gausig", "Signal Distribution1", m, rmass, width );
  RooExponential expbg( "expbg", "Background Combinatorial", m, bgpar );
  RooGaussian gausbg( "gausbg", "Background Partial Reco", m, xmass, xwidt );
  // combine them using addition
  RooAddPdf bgsum( "bgsum", "Background", RooArgList( expbg, gausbg ), fracb );
  // read values of ’m’ from a text file
  RooDataSet* data = RooDataSet::read( list, m );
  int nevr = data->numEntries();
  float nevs = nevr * fsig;
  float nevb = nevr * ( 1.0 - fsig );
  RooRealVar nsig( "nsig", "nsignal", nevs, nevs / 10.0, nevs * 10.0 );
  RooExtendPdf exsig( "exsig", "esig", gausig, nsig );
  RooRealVar nbkg( "nbkg", "nbackgr", nevb, nevb / 10.0, nevb * 10.0 );
  RooExtendPdf exbkg( "exbkg", "ebkg", bgsum, nbkg );
  RooAddPdf model( "model", "Signal + Background", RooArgList( exsig, exbkg ) );
//  RooAddPdf model( "model", "Signal + Background", RooArgList( gausig, bgsum ), fracs );
  // fit the data to the model with an unbinned maximum likelihood method
  model.fitTo( *data, Extended( kTRUE ) );
//  model.fitTo( *data );
  // print results
  cout << nfit << " results" << endl;
/*
  cout << "rmass: " << rmass.getVal() << " +/- " << rmass.getError() << endl;
  cout << "width: " << width.getVal() << " +/- " << width.getError() << endl;
  cout << "bgpar: " << bgpar.getVal() << " +/- " << bgpar.getError() << endl;
  cout << "xmass: " << xmass.getVal() << " +/- " << xmass.getError() << endl;
  cout << "xwidt: " << xwidt.getVal() << " +/- " << xwidt.getError() << endl;
  cout << "fracb: " << fracb.getVal() << " +/- " << fracb.getError() << endl;
//  cout << "fracs: " << fracs.getVal() << " +/- " << fracs.getError() << endl;
  cout << "nsig:  " << nsig .getVal() << " +/- " << nsig .getError() << endl;
  cout << "nbkg:  " << nbkg .getVal() << " +/- " << nbkg .getError() << endl;
  cout << "norm:  " << nsig .getVal() / lumi << " +/- "
                    << nsig .getError() / lumi << endl;
*/
  float sigma = width.getVal();
  ofstream os( ( nfit + ".tex" ).c_str() );
  print( rmass, os );
  print( width, os );
  print( bgpar, os );
  print( xmass, os );
  print( xwidt, os );
  print( fracb, os );
  print( nsig , os );
  print( nbkg , os );
  print( nsig , os, 1.0 / lumi );
  cout << " sigma = " << sigma << endl;
  os << "$\\sigma = " << sigma << '$' << endl;
//  int nevt = lround( data->sumEntries() );
//  fsig = fracs.getVal();
//  float esig = fracs.getError();
//  string nfit = name + dset;
//  cout << nfit << " : " << nevt * fsig << " +/- " << nevt * esig
//       << " ; sigma = " << sigma << endl;
//  cout << nfit << " : " << nevt * fsig / lumi << " +/- " << nevt * esig / lumi
//       << endl;
  // plot function
  TCanvas* can = new TCanvas( nfit.c_str(), nfit.c_str(), 800, 600 );
  can->SetLeftMargin( 0.15 );
  can->cd();
  m.setBins( nbin );
  RooPlot* mframe = m.frame();
  mframe->SetTitle( ptit.c_str() );
  mframe->SetXTitle( atit.c_str() );
  data->plotOn( mframe );
  if ( psig ) model.plotOn( mframe, Components( "gsum2" ),
                            LineStyle( kDashed ), LineColor( csig ) );
  if ( pbkg ) model.plotOn( mframe, Components( "expbg" ),
                            LineStyle( kDashed ), LineColor( cbkg ) );
  if ( pbkg ) model.plotOn( mframe, Components( "gausbg" ),
                            LineStyle( kDashed ), LineColor( cmis ) );
  model.plotOn( mframe );
  mframe->Draw();
/*
  char lstr[100];
  sprintf( lstr, "%5.2f", lumi );
  TLatex t;
  string s = "#font[";
  s += to_string( mframe->GetYaxis()->GetLabelFont() ) + "]{L}=";
  s += lstr;
  s += " fb^{-1} (#sqrt{s}=13 TeV, ";
  s += to_string( year ) + ")";
  cout << mframe->GetYaxis()->GetLabelFont() << endl;
  cout << s << endl;
  t.SetTextSize( mframe->GetYaxis()->GetLabelSize() );
  t.SetTextFont( mframe->GetYaxis()->GetLabelFont() );
  t.DrawLatex( ( xmin * 0.45 ) + ( xmax * 0.55 ), mframe->GetMaximum() * 1.01, s.c_str() );
*/
  printTitleLin( mframe, year, lumi, xmin, xmax );
  printInfo( mframe, lv, xmin, xmax );
//  float xfac = 0.65;
//  float yfac = 0.78 - ( 0.06 * lv.size() );
  yfac -= ( 0.07 + ( 0.06 * lv.size() ) );
  printLine( mframe, xmin, xmax, xfac, yfac, "Total Fit" );
  if ( psig ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Signal"           , kDashed, csig );
  if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         "Combinatorial bkg", kDashed, cbkg );
  if ( pbkg ) printLine( mframe, xmin, xmax, xfac, yfac -= 0.04,
                         misrl, kDashed, cmis );
//  nfit += ".pdf";
  can->Print( ( nfit + ".pdf" ).c_str() );
  can->Print( ( nfit + ".png" ).c_str() );

}

void fitDataset( int year, float clum, float mlum, const char* name,
                 const string& channels ) {

  set<string> chanLabels = parseLabels( channels );
  SETLOGVAL(fitijp,chanLabels);
  SETLOGVAL(fitbjp,chanLabels);
  SETLOGVAL(fitphi,chanLabels);
  SETLOGVAL(fitip2,chanLabels);
  SETLOGVAL(fitbp2,chanLabels);
  SETLOGVAL(fitups,chanLabels);
  SETLOGVAL(fitbui,chanLabels);
  SETLOGVAL(fitbud,chanLabels);
  SETLOGVAL(fitbsd,chanLabels);
  SETLOGVAL(fitbdd,chanLabels);
  SETLOGVAL(fitb0i,chanLabels);
  SETLOGVAL(fitlbi,chanLabels);
  SETLOGVAL(fitlbd,chanLabels);
  SETLOGVAL(fitbci,chanLabels);
  SETLOGVAL(fitbcd,chanLabels);
  SETLOGVAL(fitx3i,chanLabels);
  SETLOGVAL(fitx3d,chanLabels);

  vector<const InfoLabel*> ijpl;
  ijpl.push_back( new InfoLabel(
                  "#font[62]{p_{T}(#mu^{+}#mu^{-}) > 25 GeV}" ) );
  if ( fitijp )
  fitCB2ChP2   ( "list", "IncJPsi_"    , "Charmonium", year, name, clum,
                 "", "#mu^{+}#mu^{-} invariant mass [GeV]", ijpl,
                 40, 2.90, 3.30, 0.65, 0.85,
                 3.10, 0.025, 0.050, 1.3, 3.0, 0.9 );
  vector<const InfoLabel*> bjpl;
  bjpl.push_back( new InfoLabel(
                  "#font[62]{p_{T}(#mu^{+}#mu^{-}) > 20 GeV}" ) );
  bjpl.push_back( new InfoLabel(
                  "|#font[62]{y(#mu^{+}#mu^{-})| < 1.25}" ) );
  if ( fitbjp )
  fitCB2ChP2   ( "list", "BarJPsi_"    , "Charmonium", year, name, clum,
                 "", "#mu^{+}#mu^{-} invariant mass [GeV]", bjpl,
                 40, 2.90, 3.30, 0.65, 0.85,
                 3.10, 0.02, 0.05, 1.3, 3.0, 0.9 ); // only A1_lrun
  vector<const InfoLabel*> phil;
  phil.push_back( new InfoLabel(
                  "#font[62]{p_{T}(#mu^{+}#mu^{-}) > 14 GeV}" ) );
  phil.push_back( new InfoLabel(
                  "|#font[62]{y(#mu^{+}#mu^{-})| < 1.25}" ) );
  if ( fitphi )
  fitCB2ChP2   ( "list", "BarPhi_"     , "MuOnia"    , year, name, mlum,
                 "", "#mu^{+}#mu^{-} invariant mass [GeV]", phil,
                 30, 0.88, 1.18, 0.65, 0.85,
                 1.02, 0.005, 0.010, 0.5, 10.0, 0.2 );
  vector<const InfoLabel*> ip2l;
  ip2l.push_back( new InfoLabel(
                  "#font[62]{p_{T}(#mu^{+}#mu^{-}) > 18 GeV}" ) );
  if ( fitip2 )
  fitCB2ChP2   ( "list", "IncPsi2_"    , "Charmonium", year, name, clum,
                 "", "#mu^{+}#mu^{-} invariant mass [GeV]", ip2l,
                 50, 3.45, 3.95, 0.65, 0.85,
                 3.70, 0.02, 0.05, 0.5, 3.0, 0.5 );
  vector<const InfoLabel*> bp2l;
  bp2l.push_back( new InfoLabel(
                  "#font[62]{p_{T}(#mu^{+}#mu^{-}) > 10 GeV}" ) );
  bp2l.push_back( new InfoLabel(
                  "|#font[62]{y(#mu^{+}#mu^{-})| < 1.25}" ) );
  if ( fitbp2 )
  fitCB2ChP2   ( "list", "BarPsi2_"    , "Charmonium", year, name, clum,
                 "", "#mu^{+}#mu^{-} invariant mass [GeV]", bp2l,
                 50, 3.45, 3.95, 0.65, 0.85,
                 3.70, 0.02, 0.05, 0.5, 3.0, 0.5 ); // only A1_lrun
  vector<const InfoLabel*> upsl;
  upsl.push_back( new InfoLabel(
                  "#font[62]{p_{T}(#mu^{+}#mu^{-}) > 12 GeV}" ) );
  upsl.push_back( new InfoLabel(
                  "|#font[62]{y(#mu^{+}#mu^{-})| < 1.4}" ) );
  if ( fitups )
  fitCB4ChP2   ( "list", "BarUpsilon123_", "MuOnia", year, name, mlum,
                 "", "#mu^{+}#mu^{-} invariant mass [GeV]", upsl,
		 115, 8.70, 11.00, 0.65, 0.85,
                 9.45, 0.03, 0.10, 0.5, 3.0,
                10.00, 0.05, 0.5, 3.0, 10.35, 0.05, 0.5, 3.0, 0.5 );
  vector<const InfoLabel*> buil;
  buil.push_back( new InfoLabel(
                  "#font[62]{B^{#pm} #rightarrow J/#psi K^{#pm}}" ) );
  buil.push_back( new InfoLabel(
                  "#font[62]{p_{T} > 27 GeV}" ) );
  if ( fitbui )
  fitGauss2ExpE2D( "list", "InclusiveBu_", "Charmonium", year, name, clum,
                   "", "M(J/#psi K^{#pm}) [GeV]", buil,
                   "#font[62]{J/#psi K + X }",
                   100, 5.0, 6.0, 0.65, 0.85, 0.65, 0.85,
                   5.27, 0.01, 0.02, -3.0, 5.10, 0.01, 0.95, 0.5, true, true );
  vector<const InfoLabel*> budl;
  budl.push_back( new InfoLabel(
                  "#font[62]{B^{#pm} #rightarrow J/#psi K^{#pm}}" ) );
  budl.push_back( new InfoLabel(
                  "#font[62]{p_{T} > 10 GeV}" ) );
  if ( fitbud )
  fitGauss2ExpE( "list", "DisplacedBu_", "Charmonium", year, name, clum,
                 "", "M(J/#psi K^{#pm}) [GeV]", budl,
                 "#font[62]{J/#psi K + X }",
                 100, 5.0, 6.0, 0.65, 0.85,
                 5.27, 0.01, 0.02, -3.0, 5.10, 0.01, 0.95, 0.5, true, true );
  vector<const InfoLabel*> bsdl;
  bsdl.push_back( new InfoLabel(
                  "#font[62]{B^{0}_{s} #rightarrow J/#psi #Phi}" ) );
  bsdl.push_back( new InfoLabel(
                  "#font[62]{p_{T} > 10 GeV}" ) );
  if ( fitbsd )
  fitGauss2Exp ( "list", "DisplacedBs_", "Charmonium", year, name, clum,
                 "", "M(J/#psi K^{+}K^{-}) [GeV]", bsdl,
                 100, 5.0, 6.0, 0.65, 0.85,
                 5.37, 0.01, 0.02, -3.0, 0.5, true, true );
  vector<const InfoLabel*> bddl;
  bddl.push_back( new InfoLabel(
                  "#font[62]{B^{0} #rightarrow J/#psi K^{*0}}" ) );
  bddl.push_back( new InfoLabel(
                  "#font[62]{p_{T} > 10 GeV}" ) );
  if ( fitbdd )
  fitGauss2Exp ( "list", "DisplacedBd_", "Charmonium", year, name, clum,
                 "", "M(J/#psi K^{#pm}#pi^{#mp}) [GeV]", bddl,
                 100, 5.0, 6.0, 0.65, 0.85,
                 5.27, 0.01, 0.02, -3.0, 0.5, true, true );
  vector<const InfoLabel*> b0il;
  b0il.push_back( new InfoLabel(
                  "#font[62]{B^{0} #rightarrow J/#psi K^{0}_{S}}" ) );
  b0il.push_back( new InfoLabel(
                  "#font[62]{p_{T} > 10 GeV}" ) );
////  if ( fitb0i )
////  fitGauss2Exp ( "list", "InclusiveB0_", "Charmonium", year, name, clum,
////                 "", "M(J/#psi K^{+}K^{-}) [GeV]", b0il,
////                 100, 5.0, 6.0, 0.65, 0.78,
////                 5.27, 0.005, 0.020, -3.0, 0.5, true, true );
//  if ( fitb0i )
//  fitGauss2ExpG( "list", "InclusiveB0_", "Charmonium", year, name, clum,
//                 "", "M(J/#psi K^{+}K^{-}) [GeV]", b0il,
//                 "J/#psi K^{0}_{S} + X",
//                 100, 5.0, 6.0, 0.65, 0.85,
//                 5.27, 0.005, 0.020, -3.0, 5.05, 0.01, 0.95, 0.5, true, true );
//  if ( fitb0i )
//  fitGaussExpG ( "list", "InclusiveB0_", "Charmonium", year, name, clum,
//                 "", "M(J/#psi K^{+}K^{-}) [GeV]", b0il,
//                 "J/#psi K^{0}_{S} + X",
//                 100, 5.0, 6.0, 0.65, 0.85,
//                 5.27, 0.015, -3.0, 5.10, 0.01, 0.95, 0.5, true, true );
  if ( fitb0i )
  fitGaussExpE ( "list", "InclusiveB0_", "Charmonium", year, name, clum,
                 "", "M(J/#psi K^{0}_{S}) [GeV]", b0il,
                 "J/#psi K^{0}_{S} + X",
                 100, 5.0, 6.0, 0.65, 0.85,
                 5.27, 0.015, -3.0, 5.10, 0.1, 0.95, 0.5, true, true );
  vector<const InfoLabel*> lbil;
  lbil.push_back( new InfoLabel(
                  "#font[62]{#Lambda_{b} #rightarrow J/#psi #Lambda^{0}}" ) );
  lbil.push_back( new InfoLabel(
                  "#font[62]{p_{T} > 27 GeV}" ) );
  if ( fitlbi )
  fitGaussExp  ( "list", "InclusiveLambdab_", "Charmonium", year, name, clum,
                 "", "M(J/#psi #Lambda^{0}) [GeV]", lbil,
                 100, 5.0, 6.0, 0.25, 0.98,
                 5.62, 0.015, -3.0, 0.05, true, true );
  vector<const InfoLabel*> lbdl;
  lbdl.push_back( new InfoLabel(
                  "#font[62]{#Lambda_{b} #rightarrow J/#psi #Lambda^{0}}" ) );
  lbdl.push_back( new InfoLabel(
                  "#font[62]{p_{T} > 10 GeV}" ) );
  if ( fitlbd )
  fitGaussExp  ( "list", "DisplacedLambdab_", "Charmonium", year, name, clum,
                 "", "M(J/#psi #Lambda^{0}) [GeV]", lbdl,
                 100, 5.0, 6.0, 0.25, 0.98,
                 5.62, 0.015, -3.0, 0.05, true, true );
  vector<const InfoLabel*> bcil;
  bcil.push_back( new InfoLabel(
                  "#font[62]{B_{c}^{#pm} #rightarrow J/#psi #pi^{#pm}}" ) );
  bcil.push_back( new InfoLabel(
                  "#font[62]{p_{T} > 27 GeV}" ) );
  if ( fitbci )
  fitGaussExp  ( "list", "InclusiveBc_", "Charmonium", year, name, clum,
                 "", "M(J/#psi #pi^{#pm}) [GeV]", bcil,
                 100, 6.0, 7.0, 0.50, 0.98,
                 6.27, 0.05, -3.0, 0.10, true, true );
  vector<const InfoLabel*> bcdl;
  bcdl.push_back( new InfoLabel(
                  "#font[62]{B_{c}^{#pm} #rightarrow J/#psi #pi^{#pm}}" ) );
  bcdl.push_back( new InfoLabel(
                  "#font[62]{p_{T} > 10 GeV}" ) );
  if ( fitbcd )
  fitGaussExp  ( "list", "DisplacedBc_", "Charmonium", year, name, clum,
                 "", "M(J/#psi #pi^{#pm}) [GeV]", bcdl,
                 100, 6.0, 7.0, 0.50, 0.98,
                 6.27, 0.05, -3.0, 0.10, true, true );


  vector<const InfoLabel*> x3il;
  x3il.push_back( new InfoLabel(
                  "#font[62]{J/#psi #pi^{+} #pi^{-}}" ) );
  x3il.push_back( new InfoLabel(
                  "#font[62]{p_{T} > 27 GeV}" ) );
//  if ( fitx3i )
//  fitCBChP2    ( "list", "InclusiveX3872_", "Charmonium", year, name, clum,
//                 "", "M(J/#psi #pi^{+} #pi^{-}) [GeV]", x3il,
//                 40, 3.6, 4.0, 0.50, 0.98,
//                 3.69, 0.01, 1.3, 3.0, 0.1, false, true, true, true );
  if ( fitx3i )
  fitCBGChP2   ( "list", "InclusiveX3872_", "Charmonium", year, name, clum,
                 "", "M(J/#psi #pi^{+} #pi^{-}) [GeV]", x3il,
                 "#psi(2S)", "X(3872)",
                 40, 3.6, 4.0, 0.50, 0.98,
                 3.69, 0.01, 1.3, 3.0,
                 3.87, 0.01, 0.1, 0.005, true, true, true, true );
  vector<const InfoLabel*> x3dl;
  x3dl.push_back( new InfoLabel(
                  "#font[62]{J/#psi #pi^{+} #pi^{-}}" ) );
  x3dl.push_back( new InfoLabel(
                  "#font[62]{p_{T} > 10 GeV}" ) );
//  if ( fitx3d )
//  fitCBChP2    ( "list", "DisplacedX3872_", "Charmonium", year, name, clum,
//                 "", "M(J/#psi #pi^{+} #pi^{-}) [GeV]", x3dl,
//                 40, 3.6, 4.0, 0.50, 0.98,
//                 3.69, 0.01, 1.3, 3.0, 0.1, false, true, true, true );
  if ( fitx3d )
  fitCBGChP2   ( "list", "DisplacedX3872_", "Charmonium", year, name, clum,
                 "", "M(J/#psi #pi^{+} #pi^{-}) [GeV]", x3dl,
                 "#psi(2S)", "X(3872)",
                 40, 3.6, 4.0, 0.50, 0.98,
                 3.69, 0.01, 1.3, 3.0,
                 3.87, 0.01, 0.1, 0.005, true, true, true, true );
  return;
}

void histoFit( const string& datasets, const string& channels ) {
  set<string> dataLabels = parseLabels( datasets );
  SETLOGVAL(fitA,dataLabels);
  SETLOGVAL(fitB,dataLabels);
  SETLOGVAL(fitC,dataLabels);
  gROOT->ProcessLine( ".L tdrstyle.C" );
/*
  float clum;
  float mlum;
  string name;
  int year = 2018;
  name = "A1_recoHlrun", clum = 5.421429; mlum = 5.421429;
  name = "A1_recoHhrun", clum = 3.003443; mlum = 3.003443;
  name = "A2_recoHhrun", clum = 0.839398; mlum = 0.839398;
  fitCB2ChP2  ( "list", "IncJPsi_"    , "Charmonium", year, name, clum,
                "", "#mu^{+}#mu^{-} invariant mass [GeV]",
                40, 2.90, 3.30, 0.65, 0.85,
                3.10, 0.025, 0.050, 1.3, 3.0, 0.9 );
*/
  if ( fitA ) fitDataset( 2018, 13.872646,  8.766920, "A", channels );
  if ( fitB ) fitDataset( 2018,  6.863278,  6.863278, "B", channels );
  if ( fitC ) fitDataset( 2018,  6.612243,  6.280281, "C", channels );
}

