boundaryDict = {"SubsonicInletLeft":(1,0,1,0,1,18,0,-1),       "SupersonicInletLeft":(1,0,2,0,2,18,0,-1),
                "SubsonicOutletLeft":(1,0,3,0,3,18,0,-1),      "SupersonicOutletLeft":(1,0,4,0,4,18,0,-1),
                "AdiabaticNoSlipWallLeft":(1,0,5,0,5,18,0,-1), "IsotermalNoSlipWallLeft":(1,0,6,0,6,18,0,-1),

                "SubsonicInletRight":(1,0,9,0,9,18,0,1),       "SupersonicInletRight":(1,0,10,0,10,18,0,1),
                "SubsonicOutletRight":(1,0,11,0,11,18,0,1),      "SupersonicOutletRight":(1,0,12,0,12,18,0,1),
                "AdiabaticNoSlipWallRight":(1,0,13,0,13,18,0,1), "IsotermalNoSlipWallRight":(1,0,14,0,14,18,0,1),

                "SubsonicInletBottom":(0,1,0,1,1,22,-1,0),       "SupersonicInletBottom":(0,1,0,2,2,22,-1,0),
                "SubsonicOutletBottom":(0,1,0,3,3,22,-1,0),      "SupersonicOutletBottom":(0,1,0,4,4,22,-1,0),
                "AdiabaticNoSlipWallBottom":(0,1,0,5,5,22,-1,0), "IsotermalNoSlipWallBottom":(0,1,0,6,6,22,-1,0),

                "SubsonicInletTop":(0,1,0,9,9,22,1,0),        "SupersonicInletTop":(0,1,0,10,10,22,1,0),
                "SubsonicOutletTop":(0,1,0,11,11,22,1,0),      "SupersonicOutletTop":(0,1,0,12,12,22,1,0),
                "AdiabaticNoSlipWallTop":(0,1,0,13,13,22,1,0), "IsotermalNoSlipWallTop":(0,1,0,14,14,22,1,0) }
convertTimeFunctionNameToFlag = {
                                 "linearRamp":1,
                                 "exponentialRamp":2,
                                 "oscillatory":3,
                                 "multiLinearRamp":4
                                 }
