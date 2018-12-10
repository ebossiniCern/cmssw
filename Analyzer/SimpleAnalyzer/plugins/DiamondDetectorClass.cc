
#include "DiamondDetectorClass.h"




DiamondDetectorClass::DiamondDetectorClass(const edm::ParameterSet& iConfig, edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondRecHit> > tokenRecHit_input, edm::EDGetTokenT< edm::DetSetVector<CTPPSDiamondLocalTrack> > tokenLocalTrack_input)
:
SPC_mapV_(Correction_Map:: create_correction_mapV()),
tokenRecHit_          ( tokenRecHit_input),
tokenLocalTrack_          ( tokenLocalTrack_input),
RecHit_mapV_ (2),
saturationV_ (2)
{

}

DiamondDetectorClass::~ DiamondDetectorClass()
{
}


void DiamondDetectorClass::ExtractData(const edm::Event& iEvent)
{
	
	
	edm::Handle< edm::DetSetVector<CTPPSDiamondRecHit> > timingRecHit;
	edm::Handle< edm::DetSetVector<CTPPSDiamondLocalTrack> > LocalTracks;
	iEvent.getByToken( tokenRecHit_, timingRecHit );
	iEvent.getByToken( tokenLocalTrack_, LocalTracks );
	
	
	
	RecHit_mapV_.clear();
	saturationV_.clear();
	Mux_map_.clear();
	Mux_validT_map_.clear();
  
    // extract reco hit
  
	for (const auto& recHits : *timingRecHit) //rechits = array of hits in one channel
	{
		const CTPPSDiamondDetId detId( recHits.detId() );
		
		// retrieve and order all events in map. 
		
		for (const auto& recHit : recHits ) //rechit
		{
			//std::cout << "Hits in channel " << detId.channel() << ", plane " << detId.plane() << ", rp " << detId.rp()<< ", station " << detId.station()<< ", arm " << detId.arm() << std::endl;
		
			Mux_map_[std::make_pair(detId.arm(),detId.plane())]++;
			
			//put in hit map ("select hit with valid leading time and TOT")
			
			if(recHit.getT()!=0.0 && recHit.getToT()!= 0.0) 
			{
				RecHit_mapV_[detId.arm()][std::make_pair(detId.plane(),detId.channel())].push_back(recHit);
				Mux_validT_map_[std::make_pair(detId.arm(),detId.plane())]++;
			}
		
		}
	}
  
	//check for saturation
  
	for (int sector_id=0; sector_id <2; sector_id++)
	{		
		for  (int ch_number=0; ch_number < CHANNELS_X_PLANE; ch_number++)
		{
			if (RecHit_mapV_[sector_id].find( std::make_pair(PLANE_2_ID,ch_number)) != RecHit_mapV_[sector_id].end() && RecHit_mapV_[sector_id].find( std::make_pair(PLANE_3_ID,ch_number)) != RecHit_mapV_[sector_id].end())
			{	 
				if (RecHit_mapV_[sector_id][std::make_pair(PLANE_2_ID,ch_number)].at(0).getToT() > 15 && RecHit_mapV_[sector_id][std::make_pair(PLANE_3_ID,ch_number)].at(0).getToT() > 15) // double saturated
					saturationV_[sector_id] = 1;	
			}
		}
	}
	 
	//extract tracks
	int nTracks=0;
	for (const auto& locTracks : *LocalTracks) //rechits = array of hits in one channel
	{
		for (const auto& locTrack : locTracks)
		{
			nTracks++; 
			std::cout << locTrack.getX0() << std::endl;
		}
	}		
  
  std::cout << "tracks found in event: "<< nTracks << std::endl;
}

double DiamondDetectorClass::GetTime_SPC(int sector, int plane, int channel, double slope_add)
{
	double Time_raw = RecHit_mapV_[sector][std::make_pair(plane, channel)].at(0).getT();
	double ToT = RecHit_mapV_[sector][std::make_pair(plane, channel)].at(0).getToT();
	
	double C0 = SPC_mapV_[sector][std::make_pair(plane,channel)].at(0);
	double C1 = SPC_mapV_[sector][std::make_pair(plane,channel)].at(1)+slope_add;
	double C2 = SPC_mapV_[sector][std::make_pair(plane,channel)].at(2);
	double C3 = SPC_mapV_[sector][std::make_pair(plane,channel)].at(3);
	
	double Time_SPC = Time_raw - (C0 + C1*ToT + C2*ToT*ToT + C3*ToT*ToT*ToT);
 
	return Time_SPC;
	
	
}

int DiamondDetectorClass::GetSpread(int sector, int plane)
{
	
	int min=13, max=-1;
	for (const auto& rechit : RecHit_mapV_[sector])
	{
		int rh_plane = rechit.first.first;
		if (rh_plane == plane)
		{
			int rh_channel = Ch_position[rechit.first.second]; 
			if (min > rh_channel) min = rh_channel;
			if (max < rh_channel) max = rh_channel;
		}
	}
	
	
	if ((max-min)>=0) return max-min+1;
	else return 0;
	
	
}








