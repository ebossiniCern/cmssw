#ifndef BOSS__DIAMONDDETECTORCLASS__LIBGUARD
#define BOSS__DIAMONDDETECTORCLASS__LIBGUARD


// system include files
#include <memory> 

// user include files
//#include "FWCore/Framework/interface/Frameworkfwd.h"
#include "FWCore/Framework/interface/one/EDAnalyzer.h"

#include "FWCore/Framework/interface/Event.h"
#include "FWCore/Framework/interface/MakerMacros.h"
#include "FWCore/Framework/interface/ESHandle.h"
#include "DataFormats/Common/interface/DetSetVector.h"
#include "FWCore/ServiceRegistry/interface/Service.h"
#include "CommonTools/UtilAlgos/interface/TFileService.h"

#include "FWCore/ParameterSet/interface/ParameterSet.h"
#include "FWCore/Utilities/interface/InputTag.h"



#include "DataFormats/CTPPSDigi/interface/CTPPSDiamondDigi.h"
#include "DataFormats/CTPPSDetId/interface/CTPPSDiamondDetId.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondRecHit.h"
#include "DataFormats/CTPPSReco/interface/CTPPSDiamondLocalTrack.h"

#include <map>




struct Correction_Map{
    static std::vector< std::map<std::pair<int,int>,std::vector<double> > > create_correction_mapV() // plane, channel
        {
          std::vector< std::map<std::pair<int,int>,std::vector<double> > > m(2);
		  m[0][std::make_pair(0,0)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(0,1)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(0,2)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(0,3)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(0,4)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(0,5)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(0,6)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(0,7)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(0,8)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(0,9)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(0,10)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(0,11)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(1,0)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(1,1)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(1,2)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(1,3)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(1,4)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(1,5)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(1,6)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(1,7)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(1,8)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(1,9)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(1,10)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(1,11)] = std::vector<double>(4,0.0);
		  m[0][std::make_pair(2,0)] = std::vector<double> {9.60,-0.350,0,0};
		  m[0][std::make_pair(2,1)] = std::vector<double> {7.75,-0.250,0,0};
		  m[0][std::make_pair(2,2)] = std::vector<double> {8.36,-0.260,0,0};
		  m[0][std::make_pair(2,3)] = std::vector<double> {9.18,-0.324,0,0};
		  m[0][std::make_pair(2,4)] = std::vector<double> {10.4,-0.350,0,0};
		  m[0][std::make_pair(2,5)] = std::vector<double> {7.50,-0.219,0,0};
		  m[0][std::make_pair(2,6)] = std::vector<double> {9.00,-0.200,0,0};
		  m[0][std::make_pair(2,7)] = std::vector<double> {9.75,-0.258,0,0};
		  m[0][std::make_pair(2,8)] = std::vector<double> {10.1,-0.300,0,0};
		  m[0][std::make_pair(2,9)] = std::vector<double> {8.22,-0.250,0,0};
		  m[0][std::make_pair(2,10)] = std::vector<double>{8.19,-0.230,0,0};
		  m[0][std::make_pair(2,11)] = std::vector<double>{9.36,-0.283,0,0};
		  m[0][std::make_pair(3,0)] = std::vector<double> {9.29,-0.416,0,0};
		  m[0][std::make_pair(3,1)] = std::vector<double> {8.55,-0.300,0,0};
		  m[0][std::make_pair(3,2)] = std::vector<double> {9.29,-0.320,0,0};
		  m[0][std::make_pair(3,3)] = std::vector<double> {8.95,-0.300,0,0};
		  m[0][std::make_pair(3,4)] = std::vector<double> {8.35,-0.350,0,0};
		  m[0][std::make_pair(3,5)] = std::vector<double> {9.55,-0.221,0,0};
		  m[0][std::make_pair(3,6)] = std::vector<double> {9.00,-0.200,0,0};
		  m[0][std::make_pair(3,7)] = std::vector<double> {8.17,-0.350,0,0};
		  m[0][std::make_pair(3,8)] = std::vector<double> {9.34,-0.225,0,0};
		  m[0][std::make_pair(3,9)] = std::vector<double> {9.24,-0.320,0,0};
		  m[0][std::make_pair(3,10)] = std::vector<double>{8.82,-0.215,0,0};
		  m[0][std::make_pair(3,11)] = std::vector<double>{8.87,-0.325,0,0};
          
		  m[1][std::make_pair(0,0)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(0,1)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(0,2)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(0,3)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(0,4)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(0,5)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(0,6)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(0,7)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(0,8)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(0,9)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(0,10)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(0,11)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(1,0)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(1,1)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(1,2)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(1,3)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(1,4)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(1,5)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(1,6)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(1,7)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(1,8)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(1,9)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(1,10)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(1,11)] = std::vector<double>(4,0.0);
		  m[1][std::make_pair(2,0)] = std::vector<double> {0,-0.3,0,0};
		  m[1][std::make_pair(2,1)] = std::vector<double> {0,-0.3,0,0};
		  m[1][std::make_pair(2,2)] = std::vector<double> {0,-0.25,0,0};
		  m[1][std::make_pair(2,3)] = std::vector<double> {0,-0.25,0,0};
		  m[1][std::make_pair(2,4)] = std::vector<double> {0,-0.3,0,0};
		  m[1][std::make_pair(2,5)] = std::vector<double> {0,-0.2,0,0};
		  m[1][std::make_pair(2,6)] = std::vector<double> {0,-0.2,0,0};
		  m[1][std::make_pair(2,7)] = std::vector<double> {0,-0.25,0,0};
		  m[1][std::make_pair(2,8)] = std::vector<double> {0,-0.3,0,0};
		  m[1][std::make_pair(2,9)] = std::vector<double> {0,-0.3,0,0};
		  m[1][std::make_pair(2,10)] = std::vector<double>{0,-0.25,0,0};
		  m[1][std::make_pair(2,11)] = std::vector<double>{0,-0.3,0,0};
		  m[1][std::make_pair(3,0)] = std::vector<double> {0,-0.3,0,0};
		  m[1][std::make_pair(3,1)] = std::vector<double> {0,-0.3,0,0};
		  m[1][std::make_pair(3,2)] = std::vector<double> {0,-0.3,0,0};
		  m[1][std::make_pair(3,3)] = std::vector<double> {0,-0.25,0,0};
		  m[1][std::make_pair(3,4)] = std::vector<double> {0,-0.35,0,0};
		  m[1][std::make_pair(3,5)] = std::vector<double> {0,-0.25,0,0};
		  m[1][std::make_pair(3,6)] = std::vector<double> {0,-0.2,0,0};
		  m[1][std::make_pair(3,7)] = std::vector<double> {0,-0.25,0,0};
		  m[1][std::make_pair(3,8)] = std::vector<double> {0,-0.25,0,0};
		  m[1][std::make_pair(3,9)] = std::vector<double> {0,-0.3,0,0};
		  m[1][std::make_pair(3,10)] = std::vector<double>{0,-0.3,0,0};
		  m[1][std::make_pair(3,11)] = std::vector<double>{0,-0.3,0,0};
          return m;
        }

};

const std::vector<int> Ch_position {0,2,4,6,8,10,11,9,7,5,3,1};

class DiamondDetectorClass  {
	
	typedef std::map< std::pair<int,int> , std::vector<CTPPSDiamondRecHit> >  RecHit_map;
	
	
   public:
      explicit DiamondDetectorClass(const edm::ParameterSet&, edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondRecHit> >, edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondLocalTrack> > );
      ~ DiamondDetectorClass();
	  
	  void ExtractData(const edm::Event&);
	  
	  inline int GetMux(int sector, int plane) 			{return Mux_map_[std::make_pair(sector,plane)];}
	  inline int GetMuxValidT(int sector, int plane) 	{return Mux_validT_map_[std::make_pair(sector,plane)];}
	  inline bool SecIsSaturated(int sector) 			{return saturationV_[sector]!=0;}
	  
	  inline bool PairActive(int sector, int planeA, int channelA, int planeB, int channelB)   			
		{ return (RecHit_mapV_[sector].find( std::make_pair(planeA,channelA)) != RecHit_mapV_[sector].end() && RecHit_mapV_[sector].find( std::make_pair(planeB,channelB)) != RecHit_mapV_[sector].end());}
		
	  inline bool PadActive(int sector, int plane, int channel)   			
		{ return (RecHit_mapV_[sector].find( std::make_pair(plane,channel)) != RecHit_mapV_[sector].end());}
		
	  inline double GetTime(int sector, int plane, int channel)
		{return RecHit_mapV_[sector][std::make_pair(plane, channel)].at(0).getT();}	
		
	  inline std::vector<CTPPSDiamondRecHit>  GetPad(int sector, int plane, int channel) {return RecHit_mapV_[sector][std::make_pair(plane, channel)];};
		
	  double GetTime_SPC(int, int , int  , double =0.0 );
	  
	  int GetSpread(int, int);


   private:
   
   	enum Sector_id_
	{
		SECTOR_45_ID,
		SECTOR_56_ID
    };
	  
	enum Plane_id_
	{
		PLANE_0_ID,
		PLANE_1_ID,
		PLANE_2_ID,
		PLANE_3_ID
		
    };
	
    static const int CHANNELS_X_PLANE=12;
    static const int PLANES_X_DETECTOR=4;
    static const int MAX_SECTOR_NUMBER=2;
	
    std::vector<std::map<std::pair<int,int>,std::vector<double> > > SPC_mapV_;

	// ---------- objects to retrieve ---------------------------
	edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondRecHit> > tokenRecHit_;
	edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondLocalTrack> > tokenLocalTrack_;
	
	// ---------- data extracted---------------------------
	
	std::vector<RecHit_map> RecHit_mapV_;
	
	// ---------- mux map ---------------------------
	
	std::map< std::pair< int , int >, int> Mux_map_;   //arm, plane
	std::map< std::pair< int , int >, int> Mux_validT_map_;  //arm, plane
	
	//  event status flags
	  
	std::vector<int> saturationV_;
};





#endif
