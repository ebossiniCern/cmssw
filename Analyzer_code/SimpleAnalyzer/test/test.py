import FWCore.ParameterSet.Config as cms
from FWCore.ParameterSet.VarParsing import VarParsing
from RecoCTPPS.TotemRPLocal.ctppsDiamondLocalReconstruction_cff import *
import string

process = cms.Process("TIMINGSTUDY")
options = VarParsing ('analysis')

process.maxEvents = cms.untracked.PSet( input = cms.untracked.int32(-1) )
process.verbosity = cms.untracked.PSet( input = cms.untracked.int32(-1) )

# minimum of logs
process.MessageLogger = cms.Service("MessageLogger",
    statistics = cms.untracked.vstring(),
    destinations = cms.untracked.vstring('cerr'),
    cerr = cms.untracked.PSet(
        threshold = cms.untracked.string('WARNING')
    )
)


options.register ('calibFile',
				  0,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Calibration input file")

options.register ('geometryFile',
				  0,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.string,
				  "Geometry input file")

options.register ('validOOT',
				  0,
				  VarParsing.multiplicity.singleton,
				  VarParsing.varType.int,
				  "valid OOT slice")
				  
				  
options.parseArguments()
				  

# import of standard configurations
process.load('Configuration.StandardSequences.Services_cff')
process.load('FWCore.MessageService.MessageLogger_cfi')
process.load('Configuration.EventContent.EventContent_cff')
process.load('Configuration.StandardSequences.FrontierConditions_GlobalTag_cff')

# load DQM framework 
#process.load("DQM.Integration.config.environment_cfi")
#process.dqmEnv.subSystemFolder = "CTPPS"
#process.dqmEnv.eventInfoFolder = "EventInfo"
#process.dqmSaver.path = ""
#process.dqmSaver.tag = "CTPPS"

# process.source = cms.Source("NewEventStreamFileReader",
#     fileNames = cms.untracked.vstring(
#     *(
#     'file:///afs/cern.ch/user/n/nminafra/Work/public/SampicCMS/run319190_ls0130_streamPhysicsTOTEM3_StorageManager.dat',
#     )
#     )
# )
process.source = cms.Source("PoolSource",
    fileNames = cms.untracked.vstring(
    *(
    #'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v2/000/319/466/00000/D4CDC51F-9486-E811-8D08-FA163E0E63D7.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v2/000/319/466/00000/B4F35FF7-9986-E811-AD21-FA163E619822.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v2/000/319/466/00000/9260F9EF-9986-E811-B6FE-FA163EF4135C.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v2/000/319/466/00000/76EF31D1-9986-E811-A28F-FA163EEB45BD.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v2/000/319/466/00000/6EB624BF-9086-E811-8DB6-FA163E306BAD.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v2/000/319/466/00000/12EE7D0D-9A86-E811-92E2-FA163E04B59.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/FC90B27A-2EB4-E811-848F-02163E010D93.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/FC8051AF-24B4-E811-872B-FA163E37342E.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/FAF2F16E-27B4-E811-A1A0-FA163EBB6456.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/FAAB5471-22B4-E811-BF90-FA163E133B56.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/FA28F741-1EB4-E811-AFAD-FA163ECBC7B6.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/F6088F6F-27B4-E811-8686-FA163EDA7FBE.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/F432C69E-34B4-E811-B0E6-FA163E66E373.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/F42E0D44-38B4-E811-8AAF-FA163EA24D19.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/F2E3E56C-2EB4-E811-A795-FA163E81D8C9.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/F28E2C56-30B4-E811-AB3C-02163E010CA9.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/F0D6E7E9-2CB4-E811-B893-FA163ECC91B5.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/EEA82D4A-25B4-E811-BC71-FA163E2F2D4B.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/EA0263B2-27B4-E811-A02E-FA163E69FCF6.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/E4D7E361-23B4-E811-B696-FA163E198F4C.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/E07A000B-2DB4-E811-95EA-02163E01A101.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/DEDA2415-23B4-E811-8FD8-FA163EEE1BC5.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/DCF443D6-16B4-E811-94F7-FA163EEECA4B.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/DAC190BF-24B4-E811-A4DA-FA163EAD4CB1.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/DA84E925-1BB4-E811-B463-02163E01A068.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/D8146AE0-1FB4-E811-A3BC-02163E017661.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/D610EEC9-1CB4-E811-9EEC-FA163E94FF24.root',
'/store/data/Run2018D/ZeroBias/AOD/PromptReco-v2/000/322/332/00000/D4E0667A-2AB4-E811-8157-FA163E3AD924.root',
#'/store/data/Run2018B/ZeroBias1/AOD/PromptReco-v2/000/318/361/00000/20499182-0C78-E811-AC77-A4BF0112BD52.root',
#'/store/data/Run2018B/ZeroBias1/AOD/PromptReco-v2/000/318/361/00000/6E5F383F-0D78-E811-A8E6-FA163ECDEF72.root',
#'/store/data/Run2018B/ZeroBias1/AOD/PromptReco-v2/000/318/361/00000/BE61E668-0B78-E811-BD3D-FA163E60E3B4.root',
#'/store/data/Run2018B/ZeroBias1/AOD/PromptReco-v2/000/318/361/00000/EE9F884F-0A78-E811-8B9E-FA163EA96792.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/AA21D02F-A784-E811-8530-A4BF01277792.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/AA07FE1A-B284-E811-B673-FA163EBC0B5E.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A84D3A74-B084-E811-862D-FA163E51AFD3.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A6FCDC7C-AD84-E811-B784-02163E01A129.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A66307F6-AC84-E811-8F61-FA163E4C942A.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A65749B5-AD84-E811-AB6F-FA163EA6A331.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A6010C62-AD84-E811-B972-FA163E373CC4.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A4202676-AD84-E811-8E86-FA163E49296B.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A259320A-A984-E811-B227-FA163E1051E6.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A22F0E55-AC84-E811-B791-FA163EA8044F.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A04C1629-AB84-E811-A323-FA163ED286F5.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/A02A9F7A-B184-E811-8DD1-FA163EFD7D19.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/9EBD908F-AB84-E811-B035-02163E019ED0.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/9EB80A05-AD84-E811-A216-FA163EE0CB15.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/9E92D3C0-AD84-E811-93BB-FA163EA1E4FD.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/9E7EAFDF-B284-E811-985A-02163E010C30.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/9C7858F6-AF84-E811-8DE6-02163E01507D.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/9C78511E-AD84-E811-BC8D-FA163E126D94.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/989BC9F0-AC84-E811-A4D8-02163E01800D.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/98963900-AB84-E811-A42E-FA163E89DF43.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/968BA940-AD84-E811-B922-FA163ECE132A.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/94CEC45E-B084-E811-96EE-FA163EB83185.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/90C9F87E-B084-E811-8883-02163E019F71.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/90A68047-AB84-E811-A7D7-FA163ECEF311.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/904287C1-A784-E811-A52C-FA163E6369F5.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8EC1E27E-AD84-E811-8A42-FA163E862C57.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8CCDB182-AD84-E811-A059-FA163EF4E283.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8CB1EC5E-AD84-E811-B45E-FA163E73F4B9.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8C8271AE-AD84-E811-A5F3-02163E00BF95.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8C6ABA46-AD84-E811-A048-FA163E7B6908.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8C004A69-B084-E811-8A4C-FA163E55E9B1.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8ADD616A-AD84-E811-B94C-FA163E52AA20.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/86736060-AE84-E811-9E72-FA163EFF9AD1.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/862AD75A-AE84-E811-8BC5-FA163E376AE7.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/8409D191-AF84-E811-8E16-FA163EFE0B4E.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/82741E0F-A984-E811-913C-FA163ECFE6A5.root',
#'/store/data/Run2018C/ZeroBias/AOD/PromptReco-v1/000/319/337/00000/820B7559-AD84-E811-BA28-FA163EF96BC3.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/00183235-C47C-E811-ADB7-FA163E06876F.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/02C21249-BF7C-E811-9D6F-FA163E370861.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/0668D8CA-C37C-E811-A1D1-FA163E399762.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/08DBB77F-C97C-E811-92DB-FA163E4A4299.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/0C70B3EF-C47C-E811-991F-FA163E12E196.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/10792C91-BB7C-E811-8FB3-FA163EC10AF6.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/2488BF83-C77C-E811-B6EE-FA163EFDDA6B.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/3696483B-C67C-E811-A7B2-FA163E26BFDA.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/3C9CB4AE-C17C-E811-9902-02163E019FE8.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/3E1F414D-C67C-E811-90D7-FA163E160041.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/425EC7AA-C77C-E811-9458-02163E019F9D.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/44111C17-C57C-E811-8E15-FA163EBD447D.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/4ABA2487-C07C-E811-84E7-02163E01A179.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/4ED1812B-BD7C-E811-A217-FA163EA5977C.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/5088CA45-B07C-E811-BEE6-FA163EB83F5F.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/5E1F4173-C17C-E811-81F2-FA163E14F11E.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/6236C7F1-BF7C-E811-A4A4-02163E01A07C.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/6294EA44-C37C-E811-A096-FA163EFE36BC.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/6A092290-B97C-E811-B848-FA163E3A303D.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/6AB60FBD-C47C-E811-AB7A-FA163EA2F3FC.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/70B5D057-C37C-E811-A904-FA163E160041.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/72454809-C27C-E811-8CEE-FA163EDA08D4.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/7681FF63-C57C-E811-B0DC-FA163E850BF4.root',
#'/store/data/Run2018B/ZeroBias/AOD/PromptReco-v2/000/318/872/00000/786B16E7-BC7C-E811-99B7-FA163EC1D4C0.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/EAE56E3A-DED2-E711-9FB5-6C3BE5B59210.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/EAE3CE2B-43D3-E711-B538-6C3BE5B5A308.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/EA5FEF21-46D4-E711-9C11-001E67E6F8CD.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/EA4DDD67-33D4-E711-9BC2-001E677923DE.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/EA213543-32D3-E711-8653-44A842CF05CC.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E8FECE05-6CD5-E711-BF57-001E6779285C.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E8EE5A80-3FD4-E711-AD95-001E67792576.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E8EC2BD0-39D4-E711-81FF-A4BF0112BE28.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E855CD03-39D3-E711-BBFA-44A842CFC97E.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E83C20E4-6CD5-E711-A67A-A4BF01125588.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E825E4E6-3BD4-E711-983D-A4BF011254B8.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E81AFCE7-D9D2-E711-96DA-6C3BE5B5C0B0.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E80D2D70-6CD5-E711-9EAB-002590A3CA16.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E6D8BF2C-46D4-E711-B8C4-001E67E6F90E.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E6CCB959-6BD5-E711-BFF7-A4BF01125B78.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E6A79B26-EED2-E711-890A-A0369FC5DF9C.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E673918B-7BD3-E711-957A-44A842CFC97E.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E663F7B8-2AD3-E711-BBA8-44A842CFC9E6.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E6632F7F-3FD4-E711-AADA-001E67E6F4A9.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E64DF7BB-DED2-E711-B96D-484D7E8DF107.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E61E36CD-DAD3-E711-8CE8-008CFA0A59E0.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E60F6CC7-45D4-E711-90DA-A4BF01125E2E.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E4F4762E-3CD4-E711-9E9A-001E677926E2.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E4DE105E-47D4-E711-B823-001E673982E6.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E4B1A644-3FD3-E711-AAEC-008CFA1C6564.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E4ACE8A0-96D3-E711-ABA4-44A842CFD667.root',
#'/store/data/Run2017B/ZeroBias/AOD/17Nov2017-v1/20000/E4A46922-1FD3-E711-AD44-008CFA56D794.root'
    )
    )
)

from Configuration.AlCa.GlobalTag import GlobalTag
process.GlobalTag = GlobalTag(process.GlobalTag, '101X_dataRun2_HLT_v7', '')



# rechits production
process.load(options.geometryFile)
process.load('RecoCTPPS.TotemRPLocal.ctppsDiamondLocalReconstruction_cff')

process.diamondTimingAnalyzer = cms.EDAnalyzer('DiamondTimingAnalyzer',
 tagDigi = cms.InputTag("ctppsDiamondRawToDigi", "TimingDiamond"),
 tagRecHit = cms.InputTag("ctppsDiamondRecHits"),
 tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracks"),
 tagLocalTrack = cms.InputTag("ctppsDiamondLocalTracks"), 
 tagCalibrationFile = cms.string(options.calibFile),
 tagValidOOT = cms.int32(options.validOOT),
 Ntracks_Lcuts = cms.vint32([-1,1,-1,1]),
 Ntracks_Ucuts = cms.vint32([-1,6,-1,6]),
)

#process.PixelMux_Filter = cms.EDFilter('PixelMux_Filter',
# tagPixelLocalTrack = cms.InputTag("ctppsPixelLocalTracks"),
# Ntracks_Lcuts = cms.vint32([-1,1,-1,1]),
# Ntracks_Ucuts = cms.vint32([-1,3,-1,3]), 
#)

#process.TimingMux_Filter = cms.EDFilter('TimingMux_Filter',
# tagTimingLocalTrack = cms.InputTag("ctppsDiamondLocalTrack"),
# TTracksMux_Lcuts = cms.vint32([-1,-1]),
# TTracksMux_Ucuts = cms.vint32([-1,-1]), 
#)

process.TFileService = cms.Service("TFileService",
     fileName = cms.string(options.outputFile)
)


process.path = cms.Path(
   ctppsDiamondLocalReconstruction * process.diamondTimingAnalyzer
)
 


process.schedule = cms.Schedule(
process.path
)
