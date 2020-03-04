#!/usr/bin/python

import math


def one(FR_1,pt1):
  return FR_1[pt1]/(1-FR_1[pt1])


def two(FR_1, FR_2, pt1, pt2):
  return -FR_1[pt1]*FR_2[pt2]/((1-FR_1[pt1])*(1-FR_2[pt2]))


def derivate_one(FR_1,pt1):
  return 1.0/((1-FR_1[pt1])**2.0)

def derivate_first(FR_1, FR_2, pt1, pt2):
  return (1.0/((1-FR_1[pt1])**2.0))*(FR_2[pt2]/(1-FR_2[pt2]))


def derivate_second(FR_1, FR_2, pt1, pt2):
  return (1.0/((1-FR_2[pt2])**2.0))*(FR_1[pt1]/(1-FR_1[pt1]))


def derivate_both(FR_1, FR_2, pt1, pt2):
  return (2.0*FR_1[pt1])/((1-FR_1[pt1])**3.0)

def propagate_one(FR_1,pt1, unc):
  return math.sqrt((derivate_one(FR_1,pt1)**2.0)*(unc**2.0))


def propagate_first(FR_1, FR_2, pt1, pt2, unc):
  return math.sqrt((derivate_first(FR_1, FR_2, pt1, pt2)**2.0)*(unc**2.0))


def propagate_second(FR_1, FR_2, pt1, pt2, unc):
  return math.sqrt((derivate_second(FR_1, FR_2, pt1, pt2)**2.0)*(unc**2.0))


def propagate_both(FR_1, FR_2, pt1, pt2, unc):
  return math.sqrt((derivate_both(FR_1,FR_2,pt1,pt2)**2.0)*(unc**2.0))



def all(decayMode_2,tauPt_2,decayMode_1,tauPt_1,FR_1,FR_2,FR_3,up_1,down_1,up_2,down_2,up_3,down_3,myFR,myPt,up,down):
  for de2 in decayMode_2:
    for pt2 in tauPt_2:
      for de1 in decayMode_1:
        for pt1 in tauPt_1:

          if (pt2==-1 or de2 == -1):
            if (de1 == 0):
              nominal= FR_1[pt1]/(1-FR_1[pt1])
              if (myFR==1 and pt1 == myPt):
#                u = up_1[myPt]/(1-up_1[myPt])
#                d = down_1[myPt]/(1-down_1[myPt])
                u = nominal+propagate_one(FR_1,myPt,up_1[myPt])
                d = nominal-propagate_one(FR_1,myPt,down_1[myPt])
	      else:
		u = nominal
		d = nominal

            if (de1 == 1):
              nominal= FR_2[pt1]/(1-FR_2[pt1])
              if (myFR==2 and pt1 == myPt):
#                u = up_2[myPt]/(1-up_2[myPt])
#                d = down_2[myPt]/(1-down_2[myPt])
                u = nominal + propagate_one(FR_2,myPt,up_2[myPt])
                d = nominal - propagate_one(FR_2,myPt,down_2[myPt])
              else:
                u = nominal
                d = nominal

            if (de1 == 2):
              nominal= FR_3[pt1]/(1-FR_3[pt1])
              if (myFR==3 and pt1 == myPt):
#                u = up_3[myPt]/(1-up_3[myPt])
#                d = down_3[myPt]/(1-down_3[myPt])
                u = nominal + propagate_one(FR_3,myPt,up_3[myPt])
                d = nominal - propagate_one(FR_3,myPt,up_3[myPt])
              else:
                u = nominal
                d = nominal

          else:
            if (de1 == 0 and de2 ==0):
              nominal= -FR_1[pt1]*FR_1[pt2]/((1-FR_1[pt1])*(1-FR_1[pt2]))
              if (myFR == 1 and pt1 == myPt and pt2 == myPt):
                d = nominal - propagate_both(FR_1,FR_1,myPt,pt2,down_1[myPt])
                u = nominal + propagate_both(FR_1,FR_1,myPt,pt2,up_1[myPt])
              elif (myFR == 1 and pt1 == myPt):
#                d = -up_1[pt1]*up_1[pt2]/((1-up_1[pt1])*(1-up_1[pt2]))
#                u = -down_1[pt1]*down_1[pt2]/((1-down_1[pt1])*(1-down_1[pt2]))
                d = nominal - propagate_first(FR_1,FR_1,myPt,pt2,down_1[myPt])
                u = nominal + propagate_first(FR_1,FR_1,myPt,pt2,up_1[myPt])
              elif (myFR == 1 and pt2 == myPt):
#                d = -up_1[pt1]*up_1[pt2]/((1-up_1[pt1])*(1-up_1[pt2]))
#                u = -down_1[pt1]*down_1[pt2]/((1-down_1[pt1])*(1-down_1[pt2]))
                d = nominal - propagate_second(FR_1,FR_1,pt1,myPt,down_1[myPt])
                u = nominal + propagate_second(FR_1,FR_1,pt1,myPt,up_1[myPt])
              else:
                u = nominal
                d = nominal

            if (de1 == 1 and de2 ==0):
              nominal= -FR_2[pt1]*FR_1[pt2]/((1-FR_2[pt1])*(1-FR_1[pt2]))
              if (myFR == 2 and pt1 == myPt):
#                d = -up_2[pt1]*up_1[pt2]/((1-up_2[pt1])*(1-up_1[pt2]))
#                u = -down_2[pt1]*down_1[pt2]/((1-down_2[pt1])*(1-down_1[pt2]))
#                 d = two(up_2,FR_1,myPt,pt2)
#                 u = two(down_2,FR_1,myPt,pt2)
                d = nominal - propagate_first(FR_2,FR_1,myPt,pt2,down_2[myPt])
                u = nominal + propagate_first(FR_2,FR_1,myPt,pt2,up_2[myPt])
              elif (myFR == 1 and pt2 == myPt):
#                 d = two(FR_2,up_1,pt1,myPt)
#                 u = two(FR_2,down_1,pt1,myPt)
                d = nominal - propagate_second(FR_2,FR_1,pt1,myPt,down_1[myPt])
                u = nominal + propagate_second(FR_2,FR_1,pt1,myPt,up_1[myPt])
              else:
		u = nominal
		d = nominal

            if (de1 == 2 and de2 ==0):
              nominal= -FR_3[pt1]*FR_1[pt2]/((1-FR_3[pt1])*(1-FR_1[pt2]))
              if (myFR == 3 and pt1 == myPt):
#                d = two(up_3,FR_1,myPt,pt2)
#                u = two(down_3,FR_1,myPt,pt2)
                d = nominal - propagate_first(FR_3,FR_1,myPt,pt2,down_3[myPt])
                u = nominal + propagate_first(FR_3,FR_1,myPt,pt2,up_3[myPt])
              elif (myFR == 1 and pt2 == myPt):
#                d = two(FR_3,up_1,pt1,myPt)
#                u = two(FR_3,down_1,pt1,myPt)
                d = nominal - propagate_second(FR_3,FR_1,pt1,myPt,down_1[myPt])
                u = nominal + propagate_second(FR_3,FR_1,pt1,myPt,up_1[myPt])
              else:
		d = nominal
		u = nominal

            if (de1 == 0 and de2 ==1):
              nominal= -FR_1[pt1]*FR_2[pt2]/((1-FR_1[pt1])*(1-FR_2[pt2]))
#              d = -up_1[pt1]*up_2[pt2]/((1-up_1[pt1])*(1-up_2[pt2]))
#              u = -down_1[pt1]*down_2[pt2]/((1-down_1[pt1])*(1-down_2[pt2]))
	      if (myFR == 1 and pt1 == myPt):
#                d = two(up_1,FR_2,myPt,pt2)
#                u = two(down_1,FR_2,myPt,pt2)
                d = nominal - propagate_first(FR_1,FR_2,myPt,pt2,down_1[myPt])
                u = nominal + propagate_first(FR_1,FR_2,myPt,pt2,up_1[myPt])
              elif (myFR == 2 and pt2 == myPt):
#                d = two(FR_1,up_2,pt1,myPt)
#                u = two(FR_1,down_2,pt1,myPt)
                d = nominal - propagate_second(FR_1,FR_2,pt1,myPt,down_2[myPt])
                u = nominal + propagate_second(FR_1,FR_2,pt1,myPt,up_2[myPt])
              else:
                d = nominal
                u = nominal

            if (de1 == 1 and de2 ==1):
              nominal= -FR_2[pt1]*FR_2[pt2]/((1-FR_2[pt1])*(1-FR_2[pt2]))
#              d = -up_2[pt1]*up_2[pt2]/((1-up_2[pt1])*(1-up_2[pt2]))
#              u = -down_2[pt1]*down_2[pt2]/((1-down_2[pt1])*(1-down_2[pt2]))
              if (myFR == 2 and pt1 == myPt and pt2 == myPt):
                d = nominal - propagate_both(FR_2,FR_2,myPt,pt2,down_2[myPt])
                u = nominal + propagate_both(FR_2,FR_2,myPt,pt2,up_2[myPt])
              elif (myFR == 2 and pt1 == myPt):
#                d = two(up_2,FR_2,myPt,pt2)
#                u = two(down_2,FR_2,myPt,pt2)
                d = nominal - propagate_first(FR_2,FR_2,myPt,pt2,down_2[myPt])
                u = nominal + propagate_first(FR_2,FR_2,myPt,pt2,up_2[myPt])
              elif (myFR == 2 and pt2 == myPt):
#                d = two(FR_2,up_2,pt1,myPt)
#                u = two(FR_2,down_2,pt1,myPt)
                d = nominal - propagate_second(FR_2,FR_2,pt1,myPt,down_2[myPt])
                u = nominal + propagate_second(FR_2,FR_2,pt1,myPt,up_2[myPt])
              else:
                d = nominal
                u = nominal

            if (de1 == 2 and de2 ==1):
              nominal= -FR_3[pt1]*FR_2[pt2]/((1-FR_3[pt1])*(1-FR_2[pt2]))
#              d = -up_3[pt1]*up_2[pt2]/((1-up_3[pt1])*(1-up_2[pt2]))
#              u = -down_3[pt1]*down_2[pt2]/((1-down_3[pt1])*(1-down_2[pt2]))
              if (myFR == 3 and pt1 == myPt):
#                d = two(up_3,FR_2,myPt,pt2)
#                u = two(down_3,FR_2,myPt,pt2)
                d = nominal - propagate_first(FR_3,FR_2,myPt,pt2,down_3[myPt])
                u = nominal + propagate_first(FR_3,FR_2,myPt,pt2,up_3[myPt])
              elif (myFR == 2 and pt2 == myPt):
#                d = two(FR_3,up_2,pt1,myPt)
#                u = two(FR_3,down_2,pt1,myPt)
                d = nominal - propagate_second(FR_3,FR_2,pt1,myPt,down_2[myPt])
                u = nominal + propagate_second(FR_3,FR_2,pt1,myPt,up_2[myPt])
              else:
                d = nominal
                u = nominal

            if (de1 == 0 and de2 ==2):
              nominal= -FR_1[pt1]*FR_3[pt2]/((1-FR_1[pt1])*(1-FR_3[pt2]))
#              d = -up_1[pt1]*up_3[pt2]/((1-up_1[pt1])*(1-up_3[pt2]))
#              u = -down_1[pt1]*down_3[pt2]/((1-down_1[pt1])*(1-down_3[pt2]))
              if (myFR == 1 and pt1 == myPt):
#                d = two(up_1,FR_3,myPt,pt2)
#                u = two(down_1,FR_3,myPt,pt2)
                d = nominal - propagate_first(FR_1,FR_3,myPt,pt2,down_1[myPt])
                u = nominal + propagate_first(FR_1,FR_3,myPt,pt2,up_1[myPt])
              elif (myFR == 3 and pt2 == myPt):
#                d = two(FR_1,up_3,pt1,myPt)
#                u = two(FR_1,down_3,pt1,myPt)
                d = nominal - propagate_second(FR_1,FR_3,pt1,myPt,down_3[myPt])
                u = nominal + propagate_second(FR_1,FR_3,pt1,myPt,up_3[myPt])
              else:
                d = nominal
                u = nominal

            if (de1 == 1 and de2 ==2):
              nominal= -FR_2[pt1]*FR_3[pt2]/((1-FR_2[pt1])*(1-FR_3[pt2]))
#              d = -up_2[pt1]*up_3[pt2]/((1-up_2[pt1])*(1-up_3[pt2]))
#              u = -down_2[pt1]*down_3[pt2]/((1-down_2[pt1])*(1-down_3[pt2]))
	      if (myFR == 2 and pt1 == myPt):
#                d = two(up_2,FR_3,myPt,pt2)
#                u = two(down_2,FR_3,myPt,pt2)
                d = nominal - propagate_first(FR_2,FR_3,myPt,pt2,down_2[myPt])
                u = nominal + propagate_first(FR_2,FR_3,myPt,pt2,up_2[myPt])
              elif (myFR == 3 and pt2 == myPt):
#                d = two(FR_2,up_3,pt1,myPt)
#                u = two(FR_2,down_3,pt1,myPt)
                d = nominal - propagate_second(FR_2,FR_3,pt1,myPt,down_3[myPt])
                u = nominal + propagate_second(FR_2,FR_3,pt1,myPt,up_3[myPt])
              else:
                d = nominal
                u = nominal

            if (de1 == 2 and de2 ==2):
              nominal= -FR_3[pt1]*FR_3[pt2]/((1-FR_3[pt1])*(1-FR_3[pt2]))
#              d = -up_3[pt1]*up_3[pt2]/((1-up_3[pt1])*(1-up_3[pt2]))
#              u = -down_3[pt1]*down_3[pt2]/((1-down_3[pt1])*(1-down_3[pt2]))
              if (myFR == 3 and pt1 == myPt and pt2 == myPt):
                d = nominal - propagate_both(FR_3,FR_3,myPt,pt2,down_3[myPt])
                u = nominal + propagate_both(FR_3,FR_3,myPt,pt2,up_3[myPt])
              elif (myFR == 3 and pt1 == myPt):
#                d = two(up_3,FR_3,myPt,pt2)
#                u = two(down_3,FR_3,myPt,pt2)
                d = nominal - propagate_first(FR_3,FR_3,myPt,pt2,down_3[myPt])
                u = nominal + propagate_first(FR_3,FR_3,myPt,pt2,up_3[myPt])
              elif (myFR == 3 and pt2 == myPt):
#                d = two(FR_3,up_3,pt1,myPt)
#                u = two(FR_3,down_3,pt1,myPt)
                d = nominal - propagate_second(FR_3,FR_3,pt1,myPt,down_3[myPt])
                u = nominal + propagate_second(FR_3,FR_3,pt1,myPt,up_3[myPt])
              else:
                d = nominal
                u = nominal

	  down.append(d)
          up.append(u)
  return

#tau n.o. 2 can be 0 and 0 when it is not anti-isoated

tauPt_1 = [0,1,2]
decayMode_1 = [0,1,2]
tauPt_2 = [-1,0,1,2]
decayMode_2 = [-1,0,1,2]


#1 prong + zero pi

#muon
#FR_1 =[  0.391325175762 , 0.459732443094 , 0.479264855385  ]
#up_1 =[  0.433301717126 , 0.565544798934 , 0.593594165396  ]
#down_1 =[  0.349348634398 , 0.353920087254 , 0.364935545373  ]

#muon w HT and ltrk
#FR_1 = [  0.379430860281 , 0.414788752794 , 0.467154413462  ]
#up_1 = [  0.429723851154 , 0.528697896441 , 0.635150346468  ]
#down_1 = [  0.329137869408 , 0.300879609148 , 0.299158480455  ]

# mu w tight

FR_1 = [  0.297854751348 , 0.26441937685 , 0.368061035872  ]
up_1 = [  0.340997221997 , 0.351048844444 , 0.513531617866  ]
down_1 = [  0.2547122807 , 0.177789909256 , 0.222590453877  ]


#electron

#FR_1 =[  0.491324096918 , 0.631416857243 , 0.680804550648  ]
#up_1 =[  0.529964832298 , 0.69371126003 , 0.731387543807  ]
#down_1 =[  0.429029694131 , 0.569122454455 , 0.630221557488  ]

#1 prong + 1 pi

#muon
#FR_2=[  0.371998518705 , 0.33242508769 , 0.302239894867  ]
#up_2=[  0.396322295844 , 0.393072283735 , 0.390436799557  ]
#down_2=[  0.347674741566 , 0.271777891646 , 0.214042990177  ]

#muon new
#FR_2 = [  0.341382205486 , 0.304123044014 , 0.26732006669  ]
#up_2 = [  0.371731698462 , 0.381495403792 , 0.369265564016  ]
#down_2 = [  0.31103271251 , 0.226750684236 , 0.165374569365  ]

#mu w tight


FR_2 = [  0.241139158607 , 0.23736064136 , 0.162664860487  ]
up_2 = [  0.266080539084 , 0.306137581514 , 0.247501573004  ]
down_2 = [  0.216197778129 , 0.168583701206 , 0.0778281479696  ]


#FR_2 =[  0.365240395069 , 0.256029963493 , 0.493132531643  ]
#up_2 =[  0.390613999539 , 0.322272728959 , 0.651516535591  ]
#down_2 =[  0.298997629604 , 0.189787198028 , 0.334748527695  ]

#3 prong + zero pi

#muon
#FR_3=[  0.231257036328 , 0.291632771492 , 0.0783485323191  ]
#up_3=[  0.246783855425 , 0.342852490332 , 0.140959426582  ]
#down_3=[  0.215730217231 , 0.240413052652 , 0.0157376380564  ]

#muon new
#FR_3 = [  0.220947057009 , 0.327665776014 , 0.101208753884  ]
#up_3 = [  0.241304235632 , 0.400023586745 , 0.192824201598  ]
#down_3 = [  0.200589878386 , 0.255307965284 , 0.00959330616957  ]

# mu w tight

FR_3 = [  0.129225030541 , 0.106845401227 , 0.072904124856  ]
up_3 = [  0.144601573027 , 0.154727126114 , 0.152064081721  ]
down_3 = [  0.113848488056 , 0.0589636763406 , -0.00625583200928  ]

#FR_3 =[  0.242478907108 , 0.305328249931 , 0.282111555338  ]
#up_3 =[  0.259231648506 , 0.356793850224 , 0.365739058166  ]
#down_3 =[  0.191013306815 , 0.253862649638 , 0.19848405251  ]

nominal = -99
up = -99
down = -99

w = []
up = []
down = []

# one anti-iso: 1FR/(1-1FR)
# two anti-iso: -1FR*2FR/((1-1FR)*(1-2FR))


for de2 in decayMode_2:

  for pt2 in tauPt_2:

    for de1 in decayMode_1:

      for pt1 in tauPt_1:

        if (pt2==-1 or de2 == -1):
          if (de1 == 0):
            nominal= FR_1[pt1]/(1-FR_1[pt1])

          if (de1 == 1):
            nominal= FR_2[pt1]/(1-FR_2[pt1])

          if (de1 == 2):
            nominal= FR_3[pt1]/(1-FR_3[pt1])

        else:
          if (de1 == 0 and de2 ==0):
            nominal= -FR_1[pt1]*FR_1[pt2]/((1-FR_1[pt1])*(1-FR_1[pt2]))

          if (de1 == 1 and de2 ==0):
            nominal= -FR_2[pt1]*FR_1[pt2]/((1-FR_2[pt1])*(1-FR_1[pt2]))

          if (de1 == 2 and de2 ==0):
            nominal= -FR_3[pt1]*FR_1[pt2]/((1-FR_3[pt1])*(1-FR_1[pt2]))

          if (de1 == 0 and de2 ==1):
            nominal= -FR_1[pt1]*FR_2[pt2]/((1-FR_1[pt1])*(1-FR_2[pt2]))

          if (de1 == 1 and de2 ==1):
            nominal= -FR_2[pt1]*FR_2[pt2]/((1-FR_2[pt1])*(1-FR_2[pt2]))

          if (de1 == 2 and de2 ==1):
            nominal= -FR_3[pt1]*FR_2[pt2]/((1-FR_3[pt1])*(1-FR_2[pt2]))

          if (de1 == 0 and de2 ==2):
            nominal= -FR_1[pt1]*FR_3[pt2]/((1-FR_1[pt1])*(1-FR_3[pt2]))

          if (de1 == 1 and de2 ==2):
            nominal= -FR_2[pt1]*FR_3[pt2]/((1-FR_2[pt1])*(1-FR_3[pt2]))

          if (de1 == 2 and de2 ==2):
            nominal= -FR_3[pt1]*FR_3[pt2]/((1-FR_3[pt1])*(1-FR_3[pt2]))

        w.append(nominal)

print "import sys"

print "def QCDInvertedNormalizationSafetyCheck(era, searchMode, optimizationMode):"
print "    validForEra = 'Run2016'"
print "    validForSearchMode = '350to3000'"
print "    validForOptMode = ''"
print "    if not era == validForEra:"
print "        raise Exception('Error: inconsistent era, normalisation factors valid for',validForEra,'but trying to use with',era)"
print "    if not searchMode == validForSearchMode:"
print "        raise Exception('Error: inconsistent search mode, normalisation factors valid for',validForSearchMode,'but trying to use with',searchMode)"
print "    if not optimizationMode == validForOptMode:"
print "        raise Exception('Error: inconsistent optimization mode, normalisation factors valid for',validForOptMode,'but trying to use with',optimizationMode)"
print "QCDNormalization = {"
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2lt20:decayMode_2lt1'   :", w[0] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2lt20:decayMode_2lt1'   :", w[1] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2lt20:decayMode_2lt1'   :", w[2] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2lt20:decayMode_2lt1'   :", w[3] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2lt1'   :", w[4] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2lt1'   :", w[5] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2lt20:decayMode_2lt1'   :", w[6] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2lt20:decayMode_2lt1'   :", w[7] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2lt20:decayMode_2lt1'   :", w[8] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq20to40:decayMode_2lt1'   :", w[9] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2lt1'   :", w[10] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2lt1'   :", w[11] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2lt1'   :", w[12] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2lt1'   :", w[13] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2lt1'   :", w[14] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq20to40:decayMode_2lt1'   :", w[15] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2lt1'   :", w[16] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2lt1'   :", w[17] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq40to60:decayMode_2lt1'   :", w[18] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2lt1'   :", w[19] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2lt1'   :", w[20] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2lt1'   :", w[21] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2lt1'   :", w[22] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2lt1'   :", w[23] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq40to60:decayMode_2lt1'   :", w[24] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2lt1'   :", w[25] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2lt1'   :", w[26] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2gt60:decayMode_2lt1'   :", w[27] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2gt60:decayMode_2lt1'   :", w[28] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2gt60:decayMode_2lt1'   :", w[29] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2gt60:decayMode_2lt1'   :", w[30] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2lt1'   :", w[31] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2lt1'   :", w[32] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2gt60:decayMode_2lt1'   :", w[33] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2gt60:decayMode_2lt1'   :", w[34] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2gt60:decayMode_2lt1'   :", w[35] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2lt20:decayMode_2eq1to2'   :", w[36] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2lt20:decayMode_2eq1to2'   :", w[37] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2lt20:decayMode_2eq1to2'   :", w[38] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq1to2'   :", w[39] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq1to2'   :", w[40] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq1to2'   :", w[41] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2lt20:decayMode_2eq1to2'   :", w[42] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2lt20:decayMode_2eq1to2'   :", w[43] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2lt20:decayMode_2eq1to2'   :", w[44] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq1to2'   :", w[45] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq1to2'   :", w[46] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq1to2'   :", w[47] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq1to2'   :", w[48] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq1to2'   :", w[49] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq1to2'   :", w[50] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq1to2'   :", w[51] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq1to2'   :", w[52] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq1to2'   :", w[53] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq1to2'   :", w[54] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq1to2'   :", w[55] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq1to2'   :", w[56] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq1to2'   :", w[57] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq1to2'   :", w[58] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq1to2'   :", w[59] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq1to2'   :", w[60] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq1to2'   :", w[61] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq1to2'   :", w[62] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2gt60:decayMode_2eq1to2'   :", w[63] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2gt60:decayMode_2eq1to2'   :", w[64] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2gt60:decayMode_2eq1to2'   :", w[65] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq1to2'   :", w[66] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq1to2'   :", w[67] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq1to2'   :", w[68] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2gt60:decayMode_2eq1to2'   :", w[69] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2gt60:decayMode_2eq1to2'   :", w[70] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2gt60:decayMode_2eq1to2'   :", w[71] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2lt20:decayMode_2eq2to3'   :", w[72] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2lt20:decayMode_2eq2to3'   :", w[73] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2lt20:decayMode_2eq2to3'   :", w[74] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq2to3'   :", w[75] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq2to3'   :", w[76] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq2to3'   :", w[77] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2lt20:decayMode_2eq2to3'   :", w[78] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2lt20:decayMode_2eq2to3'   :", w[79] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2lt20:decayMode_2eq2to3'   :", w[80] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq2to3'   :", w[81] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq2to3'   :", w[82] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq2to3'   :", w[83] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq2to3'   :", w[84] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq2to3'   :", w[85] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq2to3'   :", w[86] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq2to3'   :", w[87] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq2to3'   :", w[88] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq2to3'   :", w[89] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq2to3'   :", w[90] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq2to3'   :", w[91] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq2to3'   :", w[92] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq2to3'   :", w[93] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq2to3'   :", w[94] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq2to3'   :", w[95] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq2to3'   :", w[96] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq2to3'   :", w[97] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq2to3'   :", w[98] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2gt60:decayMode_2eq2to3'   :", w[99] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2gt60:decayMode_2eq2to3'   :", w[100] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2gt60:decayMode_2eq2to3'   :", w[101] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq2to3'   :", w[102] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq2to3'   :", w[103] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq2to3'   :", w[104] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2gt60:decayMode_2eq2to3'   :", w[105] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2gt60:decayMode_2eq2to3'   :", w[106] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2gt60:decayMode_2eq2to3'   :", w[107] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2lt20:decayMode_2gt3'   :", w[108] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2lt20:decayMode_2gt3'   :", w[109] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2lt20:decayMode_2gt3'   :", w[110] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2lt20:decayMode_2gt3'   :", w[111] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2gt3'   :", w[112] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2gt3'   :", w[113] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2lt20:decayMode_2gt3'   :", w[114] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2lt20:decayMode_2gt3'   :", w[115] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2lt20:decayMode_2gt3'   :", w[116] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq20to40:decayMode_2gt3'   :", w[117] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2gt3'   :", w[118] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2gt3'   :", w[119] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2gt3'   :", w[120] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2gt3'   :", w[121] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2gt3'   :", w[122] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq20to40:decayMode_2gt3'   :", w[123] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2gt3'   :", w[124] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2gt3'   :", w[125] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq40to60:decayMode_2gt3'   :", w[126] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2gt3'   :", w[127] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2gt3'   :", w[128] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2gt3'   :", w[129] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2gt3'   :", w[130] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2gt3'   :", w[131] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq40to60:decayMode_2gt3'   :", w[132] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2gt3'   :", w[133] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2gt3'   :", w[134] , ","
print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2gt60:decayMode_2gt3'   :", w[135] , ","
print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2gt60:decayMode_2gt3'   :", w[136] , ","
print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2gt60:decayMode_2gt3'   :", w[137] , ","
print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2gt60:decayMode_2gt3'   :", w[138] , ","
print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2gt3'   :", w[139] , ","
print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2gt3'   :", w[140] , ","
print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2gt60:decayMode_2gt3'   :", w[141] , ","
print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2gt60:decayMode_2gt3'   :", w[142] , ","
print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2gt60:decayMode_2gt3'   :", w[143] , ","
print "    'Inclusive': 1,"
print "}"

u_1 = []
u_2 = []
u_3 = []

d_1 = []
d_2 = []
d_3 = []

for i in range(0,3):
  u_1.append(up_1[i]-FR_1[i])
  u_2.append(up_2[i]-FR_2[i])
  u_3.append(up_3[i]-FR_3[i])

  d_1.append(FR_1[i]-down_1[i])
  d_2.append(FR_2[i]-down_2[i])
  d_3.append(FR_3[i]-down_3[i])

for i in range(1,10):

  up = []
  down = []

  if (i==1):
    all(decayMode_2,tauPt_2,decayMode_1,tauPt_1,FR_1,FR_2,FR_3,u_1,d_1,u_2,d_2,u_3,d_3,1,0,up,down)
  elif (i==2):
    all(decayMode_2,tauPt_2,decayMode_1,tauPt_1,FR_1,FR_2,FR_3,u_1,d_1,u_2,d_2,u_3,d_3,1,1,up,down)
  elif (i==3):
    all(decayMode_2,tauPt_2,decayMode_1,tauPt_1,FR_1,FR_2,FR_3,u_1,d_1,u_2,d_2,u_3,d_3,1,2,up,down)
  elif (i==4):
    all(decayMode_2,tauPt_2,decayMode_1,tauPt_1,FR_1,FR_2,FR_3,u_1,d_1,u_2,d_2,u_3,d_3,2,0,up,down)
  elif (i==5):
    all(decayMode_2,tauPt_2,decayMode_1,tauPt_1,FR_1,FR_2,FR_3,u_1,d_1,u_2,d_2,u_3,d_3,2,1,up,down)
  elif (i==6):
    all(decayMode_2,tauPt_2,decayMode_1,tauPt_1,FR_1,FR_2,FR_3,u_1,d_1,u_2,d_2,u_3,d_3,2,2,up,down)
  elif (i==7):
    all(decayMode_2,tauPt_2,decayMode_1,tauPt_1,FR_1,FR_2,FR_3,u_1,d_1,u_2,d_2,u_3,d_3,3,0,up,down)
  elif (i==8):
    all(decayMode_2,tauPt_2,decayMode_1,tauPt_1,FR_1,FR_2,FR_3,u_1,d_1,u_2,d_2,u_3,d_3,3,1,up,down)
  elif (i==9):
    all(decayMode_2,tauPt_2,decayMode_1,tauPt_1,FR_1,FR_2,FR_3,u_1,d_1,u_2,d_2,u_3,d_3,3,2,up,down)


  print "EWKFakeTausNormalization = {"
  print "    'Inclusive': 1,"
  print "}"
  print "QCDPlusEWKFakeTausNormalization = {"
  print "    'Inclusive': 1,"
  print "}"
  print "QCDPlusEWKFakeTausNormalizationSystFakeWeightingVarDown_"+str(i)+" = {"
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2lt20:decayMode_2lt1'   :", down[0] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2lt20:decayMode_2lt1'   :", down[1] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2lt20:decayMode_2lt1'   :", down[2] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2lt20:decayMode_2lt1'   :", down[3] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2lt1'   :", down[4] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2lt1'   :", down[5] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2lt20:decayMode_2lt1'   :", down[6] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2lt20:decayMode_2lt1'   :", down[7] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2lt20:decayMode_2lt1'   :", down[8] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq20to40:decayMode_2lt1'   :", down[9] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2lt1'   :", down[10] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2lt1'   :", down[11] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2lt1'   :", down[12] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2lt1'   :", down[13] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2lt1'   :", down[14] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq20to40:decayMode_2lt1'   :", down[15] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2lt1'   :", down[16] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2lt1'   :", down[17] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq40to60:decayMode_2lt1'   :", down[18] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2lt1'   :", down[19] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2lt1'   :", down[20] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2lt1'   :", down[21] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2lt1'   :", down[22] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2lt1'   :", down[23] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq40to60:decayMode_2lt1'   :", down[24] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2lt1'   :", down[25] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2lt1'   :", down[26] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2gt60:decayMode_2lt1'   :", down[27] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2gt60:decayMode_2lt1'   :", down[28] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2gt60:decayMode_2lt1'   :", down[29] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2gt60:decayMode_2lt1'   :", down[30] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2lt1'   :", down[31] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2lt1'   :", down[32] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2gt60:decayMode_2lt1'   :", down[33] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2gt60:decayMode_2lt1'   :", down[34] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2gt60:decayMode_2lt1'   :", down[35] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2lt20:decayMode_2eq1to2'   :", down[36] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2lt20:decayMode_2eq1to2'   :", down[37] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2lt20:decayMode_2eq1to2'   :", down[38] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq1to2'   :", down[39] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq1to2'   :", down[40] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq1to2'   :", down[41] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2lt20:decayMode_2eq1to2'   :", down[42] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2lt20:decayMode_2eq1to2'   :", down[43] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2lt20:decayMode_2eq1to2'   :", down[44] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq1to2'   :", down[45] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq1to2'   :", down[46] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq1to2'   :", down[47] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq1to2'   :", down[48] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq1to2'   :", down[49] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq1to2'   :", down[50] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq1to2'   :", down[51] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq1to2'   :", down[52] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq1to2'   :", down[53] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq1to2'   :", down[54] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq1to2'   :", down[55] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq1to2'   :", down[56] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq1to2'   :", down[57] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq1to2'   :", down[58] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq1to2'   :", down[59] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq1to2'   :", down[60] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq1to2'   :", down[61] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq1to2'   :", down[62] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2gt60:decayMode_2eq1to2'   :", down[63] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2gt60:decayMode_2eq1to2'   :", down[64] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2gt60:decayMode_2eq1to2'   :", down[65] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq1to2'   :", down[66] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq1to2'   :", down[67] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq1to2'   :", down[68] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2gt60:decayMode_2eq1to2'   :", down[69] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2gt60:decayMode_2eq1to2'   :", down[70] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2gt60:decayMode_2eq1to2'   :", down[71] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2lt20:decayMode_2eq2to3'   :", down[72] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2lt20:decayMode_2eq2to3'   :", down[73] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2lt20:decayMode_2eq2to3'   :", down[74] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq2to3'   :", down[75] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq2to3'   :", down[76] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq2to3'   :", down[77] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2lt20:decayMode_2eq2to3'   :", down[78] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2lt20:decayMode_2eq2to3'   :", down[79] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2lt20:decayMode_2eq2to3'   :", down[80] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq2to3'   :", down[81] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq2to3'   :", down[82] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq2to3'   :", down[83] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq2to3'   :", down[84] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq2to3'   :", down[85] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq2to3'   :", down[86] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq2to3'   :", down[87] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq2to3'   :", down[88] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq2to3'   :", down[89] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq2to3'   :", down[90] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq2to3'   :", down[91] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq2to3'   :", down[92] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq2to3'   :", down[93] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq2to3'   :", down[94] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq2to3'   :", down[95] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq2to3'   :", down[96] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq2to3'   :", down[97] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq2to3'   :", down[98] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2gt60:decayMode_2eq2to3'   :", down[99] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2gt60:decayMode_2eq2to3'   :", down[100] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2gt60:decayMode_2eq2to3'   :", down[101] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq2to3'   :", down[102] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq2to3'   :", down[103] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq2to3'   :", down[104] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2gt60:decayMode_2eq2to3'   :", down[105] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2gt60:decayMode_2eq2to3'   :", down[106] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2gt60:decayMode_2eq2to3'   :", down[107] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2lt20:decayMode_2gt3'   :", down[108] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2lt20:decayMode_2gt3'   :", down[109] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2lt20:decayMode_2gt3'   :", down[110] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2lt20:decayMode_2gt3'   :", down[111] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2gt3'   :", down[112] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2gt3'   :", down[113] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2lt20:decayMode_2gt3'   :", down[114] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2lt20:decayMode_2gt3'   :", down[115] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2lt20:decayMode_2gt3'   :", down[116] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq20to40:decayMode_2gt3'   :", down[117] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2gt3'   :", down[118] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2gt3'   :", down[119] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2gt3'   :", down[120] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2gt3'   :", down[121] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2gt3'   :", down[122] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq20to40:decayMode_2gt3'   :", down[123] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2gt3'   :", down[124] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2gt3'   :", down[125] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq40to60:decayMode_2gt3'   :", down[126] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2gt3'   :", down[127] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2gt3'   :", down[128] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2gt3'   :", down[129] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2gt3'   :", down[130] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2gt3'   :", down[131] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq40to60:decayMode_2gt3'   :", down[132] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2gt3'   :", down[133] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2gt3'   :", down[134] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2gt60:decayMode_2gt3'   :", down[135] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2gt60:decayMode_2gt3'   :", down[136] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2gt60:decayMode_2gt3'   :", down[137] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2gt60:decayMode_2gt3'   :", down[138] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2gt3'   :", down[139] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2gt3'   :", down[140] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2gt60:decayMode_2gt3'   :", down[141] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2gt60:decayMode_2gt3'   :", down[142] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2gt60:decayMode_2gt3'   :", down[143] , ","
  print "    'Inclusive': 1,"
  print "}"
  print "QCDPlusEWKFakeTausNormalizationSystFakeWeightingVarUp_"+str(i)+" = {"
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2lt20:decayMode_2lt1'   :", up[0] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2lt20:decayMode_2lt1'   :", up[1] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2lt20:decayMode_2lt1'   :", up[2] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2lt20:decayMode_2lt1'   :", up[3] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2lt1'   :", up[4] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2lt1'   :", up[5] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2lt20:decayMode_2lt1'   :", up[6] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2lt20:decayMode_2lt1'   :", up[7] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2lt20:decayMode_2lt1'   :", up[8] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq20to40:decayMode_2lt1'   :", up[9] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2lt1'   :", up[10] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2lt1'   :", up[11] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2lt1'   :", up[12] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2lt1'   :", up[13] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2lt1'   :", up[14] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq20to40:decayMode_2lt1'   :", up[15] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2lt1'   :", up[16] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2lt1'   :", up[17] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq40to60:decayMode_2lt1'   :", up[18] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2lt1'   :", up[19] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2lt1'   :", up[20] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2lt1'   :", up[21] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2lt1'   :", up[22] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2lt1'   :", up[23] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq40to60:decayMode_2lt1'   :", up[24] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2lt1'   :", up[25] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2lt1'   :", up[26] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2gt60:decayMode_2lt1'   :", up[27] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2gt60:decayMode_2lt1'   :", up[28] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2gt60:decayMode_2lt1'   :", up[29] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2gt60:decayMode_2lt1'   :", up[30] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2lt1'   :", up[31] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2lt1'   :", up[32] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2gt60:decayMode_2lt1'   :", up[33] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2gt60:decayMode_2lt1'   :", up[34] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2gt60:decayMode_2lt1'   :", up[35] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2lt20:decayMode_2eq1to2'   :", up[36] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2lt20:decayMode_2eq1to2'   :", up[37] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2lt20:decayMode_2eq1to2'   :", up[38] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq1to2'   :", up[39] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq1to2'   :", up[40] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq1to2'   :", up[41] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2lt20:decayMode_2eq1to2'   :", up[42] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2lt20:decayMode_2eq1to2'   :", up[43] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2lt20:decayMode_2eq1to2'   :", up[44] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq1to2'   :", up[45] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq1to2'   :", up[46] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq1to2'   :", up[47] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq1to2'   :", up[48] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq1to2'   :", up[49] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq1to2'   :", up[50] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq1to2'   :", up[51] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq1to2'   :", up[52] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq1to2'   :", up[53] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq1to2'   :", up[54] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq1to2'   :", up[55] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq1to2'   :", up[56] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq1to2'   :", up[57] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq1to2'   :", up[58] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq1to2'   :", up[59] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq1to2'   :", up[60] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq1to2'   :", up[61] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq1to2'   :", up[62] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2gt60:decayMode_2eq1to2'   :", up[63] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2gt60:decayMode_2eq1to2'   :", up[64] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2gt60:decayMode_2eq1to2'   :", up[65] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq1to2'   :", up[66] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq1to2'   :", up[67] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq1to2'   :", up[68] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2gt60:decayMode_2eq1to2'   :", up[69] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2gt60:decayMode_2eq1to2'   :", up[70] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2gt60:decayMode_2eq1to2'   :", up[71] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2lt20:decayMode_2eq2to3'   :", up[72] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2lt20:decayMode_2eq2to3'   :", up[73] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2lt20:decayMode_2eq2to3'   :", up[74] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq2to3'   :", up[75] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq2to3'   :", up[76] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2eq2to3'   :", up[77] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2lt20:decayMode_2eq2to3'   :", up[78] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2lt20:decayMode_2eq2to3'   :", up[79] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2lt20:decayMode_2eq2to3'   :", up[80] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq2to3'   :", up[81] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq2to3'   :", up[82] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2eq2to3'   :", up[83] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq2to3'   :", up[84] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq2to3'   :", up[85] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2eq2to3'   :", up[86] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq2to3'   :", up[87] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq2to3'   :", up[88] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2eq2to3'   :", up[89] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq2to3'   :", up[90] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq2to3'   :", up[91] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2eq2to3'   :", up[92] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq2to3'   :", up[93] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq2to3'   :", up[94] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2eq2to3'   :", up[95] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq2to3'   :", up[96] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq2to3'   :", up[97] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2eq2to3'   :", up[98] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2gt60:decayMode_2eq2to3'   :", up[99] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2gt60:decayMode_2eq2to3'   :", up[100] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2gt60:decayMode_2eq2to3'   :", up[101] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq2to3'   :", up[102] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq2to3'   :", up[103] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2eq2to3'   :", up[104] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2gt60:decayMode_2eq2to3'   :", up[105] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2gt60:decayMode_2eq2to3'   :", up[106] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2gt60:decayMode_2eq2to3'   :", up[107] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2lt20:decayMode_2gt3'   :", up[108] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2lt20:decayMode_2gt3'   :", up[109] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2lt20:decayMode_2gt3'   :", up[110] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2lt20:decayMode_2gt3'   :", up[111] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2gt3'   :", up[112] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2lt20:decayMode_2gt3'   :", up[113] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2lt20:decayMode_2gt3'   :", up[114] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2lt20:decayMode_2gt3'   :", up[115] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2lt20:decayMode_2gt3'   :", up[116] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq20to40:decayMode_2gt3'   :", up[117] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2gt3'   :", up[118] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq20to40:decayMode_2gt3'   :", up[119] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2gt3'   :", up[120] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2gt3'   :", up[121] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq20to40:decayMode_2gt3'   :", up[122] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq20to40:decayMode_2gt3'   :", up[123] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2gt3'   :", up[124] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq20to40:decayMode_2gt3'   :", up[125] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2eq40to60:decayMode_2gt3'   :", up[126] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2gt3'   :", up[127] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2eq40to60:decayMode_2gt3'   :", up[128] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2gt3'   :", up[129] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2gt3'   :", up[130] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2eq40to60:decayMode_2gt3'   :", up[131] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2eq40to60:decayMode_2gt3'   :", up[132] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2gt3'   :", up[133] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2eq40to60:decayMode_2gt3'   :", up[134] , ","
  print "    'tauPt_1lt40:decayMode_1lt2:tauPt_2gt60:decayMode_2gt3'   :", up[135] , ","
  print "    'tauPt_1eq40to60:decayMode_1lt2:tauPt_2gt60:decayMode_2gt3'   :", up[136] , ","
  print "    'tauPt_1gt60:decayMode_1lt2:tauPt_2gt60:decayMode_2gt3'   :", up[137] , ","
  print "    'tauPt_1lt40:decayMode_1eq2to3:tauPt_2gt60:decayMode_2gt3'   :", up[138] , ","
  print "    'tauPt_1eq40to60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2gt3'   :", up[139] , ","
  print "    'tauPt_1gt60:decayMode_1eq2to3:tauPt_2gt60:decayMode_2gt3'   :", up[140] , ","
  print "    'tauPt_1lt40:decayMode_1gt3:tauPt_2gt60:decayMode_2gt3'   :", up[141] , ","
  print "    'tauPt_1eq40to60:decayMode_1gt3:tauPt_2gt60:decayMode_2gt3'   :", up[142] , ","
  print "    'tauPt_1gt60:decayMode_1gt3:tauPt_2gt60:decayMode_2gt3'   :", up[143] , ","
  print "    'Inclusive': 1,"
  print "}"
