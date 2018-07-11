float yround( float x, int s, bool c ) {
  int n = lround( x );
  int l = 0;
  int m = pow( 10, s );
  while ( x > m ) {
    x /= 10;
    ++l;
  }
  int d = x / 10;
  int a = d * 10;
  int u = x - a;
  if ( c ) {
    if ( u > 5 ) u = 10;
    else
    if ( u > 2 ) u = 5;
    else
    if ( u > 0 ) u = 2;
  }
//  cout << ( a + u ) * pow( 10, l ) << endl;
  return ( a + u ) * pow( 10, l );
}

void hPlot( TDirectory* d, const char* n ) {
  TH1F* h = dynamic_cast<TH1F*>( d->Get( n ) );
  float m = 0.0;
  int b = h->GetNbinsX();
  int i;
  for ( i = 1; i < b; ++i ) {
    int y = h->GetBinContent( i );
    if ( y > m ) m = y;
  }
  h->GetYaxis()->SetRangeUser( 0.0, yround( m * 1.05, 2, true ) );
  TCanvas* c = new TCanvas( n, n, 800, 600 );
  c->cd();
  h->Draw();
  string f = string( n ) + ".pdf";
  c->Print( f.c_str() );
  return;
}

void hPlog( TDirectory* d, const char* n ) {
  TH1F* h = dynamic_cast<TH1F*>( d->Get( n ) );
  float u = 0.0;
  float l = 1.0e+20;
  int b = h->GetNbinsX();
  int i;
  for ( i = 1; i < b; ++i ) {
    int y = h->GetBinContent( i );
    if ( y > u ) u = y;
    if ( y < l ) l = y;
  }
  u *= 3.0;
  l /= 3.0;
  if ( l < 1.0 ) l = 1.0;
  h->GetYaxis()->SetRangeUser( yround( l, 1, false ), yround( u, 1, false ) );
  TCanvas* c = new TCanvas( n, n, 800, 600 );
  c->cd();
  h->Draw();
  c->SetLogy();
  string f = string( n ) + ".pdf";
  c->Print( f.c_str() );
  return;
}

void plot() {

  gStyle->SetOptStat ( false );
  gStyle->SetOptTitle( false );

  static TDirectory* d =  dynamic_cast<TDirectory*>( _file0->Get( "bphHistoSpecificDecay" ) );
  static TDirectory* e =  dynamic_cast<TDirectory*>( _file1->Get( "bphHistoSpecificDecay" ) );

  hPlot( e, "massTBPhi"    );
//  hPlot( d, "massTBJPsi"   );
//  hPlot( d, "massTBPsi2"   );
//  hPlot( e, "massTBUps123" );

//  hPlot( d, "massTDBu"     );
//  hPlot( d, "massTDBd"     );
//  hPlot( d, "massTDBs"     );
//  hPlog( d, "ctauTIBu"     );

//  hPlot( d, "mfitTDBu"     );
//  hPlot( d, "mfitTDBd"     );
//  hPlot( d, "mfitTDBs"     );
//  hPlog( d, "ctauTBJPsi"   );

}
