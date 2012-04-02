scattering(){


string normal="Daten/testScatteringOff.txt";
string scat="Daten/testScatteringOn_2.txt";

string treename_normal="normal1e+07";
string treename_scat="scat1e+07";

TFile *file = new TFile("Daten/Scattering/testScattering.root","RECREATE");

ifstream Normal;
ifstream Scat;

Normal.open(normal.c_str());
Scat.open(scat.c_str());

double energy_normal;
double final_energy_normal;
double r_normal;
double x_normal;
double y_normal;
double z_normal;
double t_normal;
double theta_normal;
double phi_normal;

double energy_scat;
double final_energy_scat;
double r_scat;
double x_scat;
double y_scat;
double z_scat;
double t_scat;
double theta_scat;
double phi_scat;

TTree *normal_tree = new TTree(treename_normal.c_str(),treename_normal.c_str()); 
TTree *scat_tree = new TTree(treename_scat.c_str(),treename_scat.c_str());

normal_tree->Branch("energy_normal",&energy_normal,"energy_normal/D");
normal_tree->Branch("final_energy_normal",&final_energy_normal,"final_energy_normal/D");
normal_tree->Branch("r_normal",&r_normal,"r_normal/D");
normal_tree->Branch("x_normal",&x_normal,"x_normal/D");
normal_tree->Branch("y_normal",&y_normal,"y_normal/D");
normal_tree->Branch("z_normal",&z_normal,"z_normal/D");
normal_tree->Branch("t_normal",&t_normal,"t_normal/D");
normal_tree->Branch("theta_normal",&theta_normal,"theta_normal/D");
normal_tree->Branch("phi_normal",&phi_normal,"phi_normal/D");

scat_tree->Branch("energy_scat",&energy_scat,"energy_scat/D");
scat_tree->Branch("final_energy_scat",&final_energy_scat,"final_energy_scat/D");
scat_tree->Branch("r_scat",&r_scat,"r_scat/D");
scat_tree->Branch("x_scat",&x_scat,"x_scat/D");
scat_tree->Branch("y_scat",&y_scat,"y_scat/D");
scat_tree->Branch("z_scat",&z_scat,"z_scat/D");
scat_tree->Branch("t_scat",&t_scat,"t_scat/D");
scat_tree->Branch("theta_scat",&theta_scat,"theta_scat/D");
scat_tree->Branch("phi_scat",&phi_scat,"phi_scat/D");

//cout<<normal.c_str()<<endl;

while(Scat.good()){

	Normal>>energy_normal>>final_energy_normal>>r_normal>>x_normal>>y_normal>>z_normal>>t_normal>>theta_normal>>phi_normal;
	Scat>>energy_scat>>final_energy_scat>>r_scat>>x_scat>>y_scat>>z_scat>>t_scat>>theta_scat>>phi_scat;

	normal_tree->Fill();
	scat_tree->Fill();

}


normal_tree->Write();
scat_tree->Write();

file->Close();

















}
