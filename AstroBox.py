# -*- coding: UTF-8 -*-
# ----------------------------------------------------------------------------
def intro():
  """
    Moduł funkcji do obliczeń w astronomii sferycznej i klasycznej,
    korzysta z modułu do obliczeń matematycznych numpy
    * Użycie: "from AstroBox import *"
    * Lista funkcji w module: dir()
    * Opis funkcji "foo" z przykładem: opis(foo), doc(foo)
    * Uruchomienie wybranych funkcji: test()
    * Uaktualnienie modułu po zmianach w powłoce Pythona:
      from importlib import reload
      reload(AstroBox)
    * Usunięcie tego opisu: dodaj komentarz do ostatniego wiersza w pliku
      # intro(doc)  
     Krzysztof Goździewski, 2009-2026, wersja 09.04.2026
  """
#----------------------------------------------------------------------------
#
# testowanie modułu: python ./astrobox.py -v
# (niektóre wyniki pozornie wykazują błędy, które wynikają z zaokrągleń)      
#
#----------------------------------------------------------------------------
import numpy as np
import inspect
# ----------------------------------------------------------------------------
# stałe fundamentalne
# ----------------------------------------------------------------------------
# zero ,,maszynowe'', dokładność 1+eps=1
eps    = 2.26e-16
# czynnik konwersji stopni na radiany
d2r    = np.pi/180.0
# data juliańska epoki fundamentalnej J2000.0, January 1.5, 2000 (TT)
JD2000 = 2451545.0
# stulecie juliańskie
cy     = 36525.0

# ----------------------------------------------------------------------------
# Konwersje katów pomiędzy różnymi miarami oraz wypisywanie kątów z mianami
# ----------------------------------------------------------------------------
def wrap(x):
    return (x + np.pi) % (2*np.pi) - np.pi

def _angle_to_dms( rad, to_deg=True ):
    """
      zamiana radianów na (stopnie/godziny, minuty, sekundy), sekundy ∈ [0,60).
      parametry:
        rad — kąt w radianach
        to_deg — True dla stopni (180=π), False dla godzin (12h=π)
      zwraca:
        (deg, min, sec) — stopnie/godziny, minuty, sekundy
        
      przykład:
      >>> _angle_to_dms(np.pi/2, True)
      (90.0, 0.0, 0.0)
      >>> _angle_to_dms(np.pi, True)
      (180.0, 0.0, 0.0)
      >>> round(_angle_to_dms(0.12345)[2], 2)
      24.63
      >>> _angle_to_dms(np.pi, False)
      (12.0, 0.0, 0.0)
    """
    factor = 180.0/np.pi if to_deg else 12.0/np.pi
    angle = rad * factor
    
    # poprawka dla wartości ujemnych: divmod(-20.5, 1) = (-21, 0.5)
    # ale chcemy: -20° 30' 0" = (-20, 30, 0)
    sign = -1 if angle < 0 else 1
    angle_abs = np.abs(angle)
    
    deg, remainder = divmod( angle_abs, 1 )
    min_, sec = divmod(remainder * 60, 1)
    
    return int(deg)*sign, int(min_), sec*60    

def hr2rad( hr, min, sec ):
   """
     zamiana kąta w mierze czasowej na łukową
     rad = hr2rad( hr, min, sec )
     
     przykład:
     >>> round( hr2rad( 20, 30, 30 ), 9 )
     5.369069111
   """
   return( hr+ min/60.0 + sec/3600.0 )*np.pi/12.0;

def deg2rad( deg, min, sec ):
   """
     zamiana kąta w mierze stopniowej na łukową
     rad = deg2rad( deg, min, sec )
     
     przykład:
     >>> round( deg2rad( 20, 30, 30 ), 9 )
     0.357937941
   """
   return ( deg + min*1.0/60.0 + sec/3600.0 )*np.pi/180.0;

def rad2deg( rad ):
   """
     zamiana kąta w mierze łukowej na stopniową
     deg, min, sec = rad2deg( rad )
     
     przykład:
     >>> rad2deg( 0.123456 )
     (7.0, 4.0, 24.627920047527176)
     >>> printdeg( 0.123456 )
     7  4' 24.63"
   """
   return _angle_to_dms(rad, to_deg=True)

def rad2hr( rad ):
   """
     zamiana kąta w mierze łukowej na czasową
     hr, min, sec = rad2hr( rad )
     
     przykład:
     >>> rad2hr( 0.123456 )
     (0.0, 28.0, 17.64186133610181)
     >>> printhr( 0.123456 )
     0h 28m  17.64s
   """
   return _angle_to_dms( rad, to_deg=False )

def day2hr( dt ):
   """
     zamiana interwału w mierze dziennej na czasową
     day, hour, min, sec = day2hr( dt )
     
     przykład:
     >>> day2hr( 1.5 )
     (1, 12, 0, 0.0)
     >>> day2hr( 2.25 )
     (2, 6, 0, 0.0)
     >>> day2hr( 0.123456 )
     (0, 2, 57, 41.856000000000004)
   """
   day    = int( dt )
   fraction = (dt - day) * 24.0
   hour   = int( fraction )
   minute = int( (fraction - hour) * 60.0 )
   sec    = (fraction - hour - minute/60.0) * 3600.0
   return ( day, hour, minute, sec )

def print1801( angle, name='', fsec=2 ):
   """
     wypisuje wartość kąta angle [rad] w mierze stopniowej
     w zakresie [-180,180] stopni
     print180( angle, name='', fsec=2 )
     
     przykład:
     >>> print180( 3*np.pi/2. )
     -90  0' 0.00"
   """
   if (angle>np.pi):
     angle = angle-2*np.pi
   deg, min, sec = _angle_to_dms(angle, to_deg=True)
   if ( len(name)>0 ):
     format = "%s = %3.0f %2.0f' %." + str(fsec) + "f\" "
     print ( format % \
          ( name, deg, np.abs(min), np.abs(sec) ) )
   else:
     format = "%3.0f %2.0f' %." + str(fsec) + "f\" "
     print ( format % ( deg, np.abs(min), np.abs(sec) ) )
     
  

def printdeg1( angle, name='', fsec=2 ):
   """
     wypisuje wartość kąta angle [rad] o nazwie name w mierze stopniowej
     parameter domyślny fsec określa ilość cyfr znaczących w wyniku
     printdeg( angle, name='', fsec=2 )
     
     przykład:
     >>> printdeg( 3*np.pi/2. )
     270  0' 0.00"
   """
   deg, min, sec = _angle_to_dms(angle, to_deg=True)
   if ( len(name)>0 ):
     format = "%s = %3.0f %2.0f' %." + str(fsec) + "f\" "
     print( format % (name, deg, min, sec))
   else:
     format = "%3.0f %2.0f' %." + str(fsec) + "f\" "
     print( format % (deg, min, sec))
     

def print180( angle, name='', fsec=2 ):
   """
     wypisuje wartość kąta angle [rad] w mierze stopniowej
     w zakresie [-180,180] stopni
     print180( angle, name='', fsec=2 )
     przykład:
     >>> print180( 3*np.pi/2. )
     -90  0' 0.00"
   """
   # redukcja kąta do zakresu [-π, π] z zachowaniem znaku
   reduced = np.mod( angle + np.pi, 2*np.pi ) - np.pi
   deg, min_, sec = _angle_to_dms( reduced, to_deg=True )
   
   if len(name) > 0:
     format = "%s = %3.0f %2.0f' %." + str(fsec) + "f\" "
     print( format % (name, deg, min_, sec) )
   else:
     format = "%3.0f %2.0f' %." + str(fsec) + "f\" "
     print(format % (deg, min_, sec))

def printdeg(angle, name='', fsec=2):
   """
     wypisuje wartość kąta angle [rad] o nazwie name w mierze stopniowej
     parameter domyślny fsec określa ilość cyfr znaczących w wyniku
     printdeg( angle, name='', fsec=2 )
     przykład:
     >>> printdeg( 3*np.pi/2. )
     270  0' 0.00"
   """
   deg, min_, sec = _angle_to_dms(angle, to_deg=True)
   if len(name) > 0:
     format = "%s = %3.0f %2.0f' %." + str(fsec) + "f\" "
     print(format % (name, deg, min_, sec))
   else:
     format = "%3.0f %2.0f' %." + str(fsec) + "f\" "
     print(format % (deg, min_, sec))     

def printhr( angle, name='', fsec=2 ):
   """
     wypisuje wartość kąta angle [rad] o nazwie name w mierze czasowej
     parameter domyślny fsec określa ilość cyfr znaczących w wyniku
     printhr( angle, name='', fsec=2 )
     
     przykład:
     >>> printhr( 3*np.pi/2. + 0.1 )
      18h 22m  55.10s
   """
   
   angle = np.fmod( angle , 2*np.pi )
    
   deg, min, sec = _angle_to_dms(angle, to_deg=False)
   if ( len(name)>0 ):
     format = "%s = %3.0fh %2.0fm %." + str(fsec) + "fs "
     print( format % ( name, deg, min, sec ) )
   else:
     format = "%3.0fh %2.0fm %." + str(fsec) + "fs "
     print( format % ( deg, min, sec ) )

def printfhr( dt, name='' ):
   """
     zamienia datę juliańską na kalendarzową z dokładnością do minut
     printfhr( dt, name='' )
     dt - data juliańska
     przykład:
     >>> printfhr( 2450000.1, name='' )
     1995 10  9 14h 24.0m
   """
   date = JDtoDate( dt )
   yr   = date[0]
   mo   = date[1]
   day, hr, min = day2hr( date[2] );
   if ( len(name)>0 ):
     print("%s = %4d %2d %2d %2.0fh %3.1fm " % ( name, yr, mo, day, hr, min ))
   else:
     print("%4d %2d %2d %2.0fh %3.1fm " % ( yr, mo, day, hr, min ) )
     
# ----------------------------------------------------------------------------
# Przemiana współrzędnych
# ----------------------------------------------------------------------------
#
#   Własności funkcji arctan2(x1,x2), gdy tan=x1/x2:
#   x1 x2 arctan2(x1,x2)
#   +/- 0 +0 +/- 0
#   +/- 0 -0 +/- pi
#   > 0 +/-inf +0 / +pi
#   < 0 +/-inf -0 / -pi
#   +/-inf +inf +/- (pi/4)
#   +/-inf -inf +/- (3*pi/4)
#
#
# ----------------------------------------------------------------------------
# przemiana współrzędnych równikowe I <-> horyzontalne
# ----------------------------------------------------------------------------
def eqI2hor( phi, t, delta ):
   """
     współrzędne równikowe I (alpha,delta) na współrzędne horyzontalne (A,h)
     A, h = eqI2hor( phi, t, delta )
     phi jest szerokością geograficzną (geodezyjną) obserwatora
     wszystkie kąty są podane w radianach
     na podstawie wzoru cosinusów dla boków trójkąta paralaktycznego
     
     przykład:
     >>> phi=deg2rad( 60,0,0 )
     >>> printdeg( phi )
      59 60' 0.00"
     >>> phi=deg2rad( 60,2,0 )
     >>> printdeg( phi )
      60  2' 0.00"
     >>> delta=deg2rad(12,5,5)
     >>> printdeg( delta )
      12  5' 5.00"
     >>> t=hr2rad(12+9.6,0,0)
     >>> printhr( t )
     21h 36m  0.00s
     >>> A, h = eqI2hor( phi, t, delta )
     >>> t, delta = hor2eqI( phi, A, h )
     >>> printhr( t )
     21h 35m  60.00s
     >>> printhr( delta )
     >>> printdeg( delta )
     12  5' 5.00"
   """
   sh = np.sin(delta)*np.sin(phi) + np.cos(delta)*np.cos(phi)*np.cos(t)
   if ( np.abs(np.abs(phi)-np.pi/2) <eps ):
      print("# <eqI2hor>: obserwator zbyt blisko bieguna!" )
      return( 0.0, 0.0 )
   if ( np.abs((np.abs(sh)-1.0))<eps ):
      print("# <eqI2hor>: wysokość zbyt bliska +/-90!" )
      return( 0.0, 0.0 )
   h  = np.arcsin( sh );
   cA = (np.sin(delta)-np.sin(h)*np.sin(phi))/(np.cos(h)*np.cos(phi));
   sA = np.sin(t)*np.cos(delta)/np.cos(h);
   A  = np.arctan2( sA,cA );
   return( A, h );

def hor2eqI( phi, A, h ):
   """
     współrzędne horyzontalne (A,h) -> współrzędne równikowe I (t,delta)
     t, delta = hor2eqI( phi, A, h )
     phi - szerokość geograficzna
     wszystkie kąty muszą być/są podane w radianach
     na podstawie wzoru cosinusów dla boków trójkąta paralaktycznego
     
     przykład:
     >>> A = deg2rad( 224, 42, 9.01 )
     >>> h = deg2rad( 35, 12, 21.04)
     >>> phi = deg2rad( 60, 2, 0)
     >>> t, delta = hor2eqI( phi, A, h )
     >>> printhr( t )
     21h 35m  60.00s
     >>> printdeg( delta )
     12  5' 5.00"
   """
   sd = np.sin(h)*np.sin(phi) + np.cos(h)*np.cos(phi)*np.cos(A)
   if ( np.abs(np.abs(phi)-np.pi/2) <eps ):
      print("# <hor2eqI>: obserwator zbyt blisko bieguna!" )
      return( 0.0, 0.0 )
   if ( np.abs((np.abs(sd)-1.0))<eps ):
      print("# <hor2eqI>: deklinacja zbyt bliska +/-90!" )
      return( 0.0, 0.0 )
   delta  = np.arcsin( sd );
   ct = (np.sin(h)-np.sin(phi)*np.sin(delta))/(np.cos(phi)*np.cos(delta));
   st = np.sin(A)*np.cos(h)/np.cos(delta);
   t  = np.arctan2( st,ct );
   return( t, delta );

# ---------------------------------------------------------------------------
# przemiana współrzędnych równikowe II <-> horyzontalne
# ---------------------------------------------------------------------------
def eqII2hor( phi, LST, alpha, delta ):
   """
     współrzędne równikowe II (alpha,delta) na współrzędne horyzontalne (A,h)
     phi jest szerokoscią geograficzną (geodezyjną) obserwatora
     A, h = eqII2hor( phi, LST, alpha, delta )
     phi - szerokość geograficzna (geodezyjna)
     LST - miejscowy czas gwiazdowy
     alpha, delta - współrzędne równikowe
     wszystkie kąty muszą być/są podane w radianach
     funkcja wykorzystuje związek czasu gwiazdowego z rektascensją i kątem
     godzinnym oraz przemianę równikowe I -> horyzontalne
     przykład:
     >>> phi=deg2rad(53,25,25)
     >>> LST=hr2rad(9,0,0)
     >>> alpha=hr2rad(5,10,10)
     >>> delta=deg2rad(23,30,0)
     >>> A,h=eqII2hor( phi, LST, alpha, delta)
     >>> printdeg(A,"azymut")
     azymut = 101 35' 20.72"
     >>> printdeg(h,"wysokość")
     wysokość =  37 53' 30.98"
   """
   t    = LST - alpha
   A, h = eqI2hor( phi, t, delta )
   return( A, h );

def hor2eqII( phi, LST, A, h ):
   """
     współrzędne horyzontalne (A,h) -> współrzędne równikowe II (alpha,delta)
     alpha, delta = hor2eqII( phi, LST, A, h )
     phi - szerokość geograficzna
     LST - miejscowy czas gwiazdowy
     A - azymut
     h - wysokość
     wszystkie kąty muszą być/są podane w radianach
     funkcja wykorzystuje związek czasu gwiazdowego z rektascensją i kątem
     godzinnym oraz przemianę równikowe I <- horyzontalne
     przykład:
     >>> phi=deg2rad(53,25,25)
     >>> LST=hr2rad(9,0,0)
     >>> A=deg2rad(101,35,20.72)
     >>> h=deg2rad(37,53,30.98)
     >>> α,δ=hor2eqII( phi, LST, A, h )
     >>> printhr( α, "α" )
     α =   5h 10m  10.00s
     >>> printdeg( δ, "δ" )
     δ =  23 29' 60.00"
   """
   t, delta = hor2eqI( phi, A, h )
   return( LST-t, delta );

# ----------------------------------------------------------------------------
# przemiana współrzędnych równikowe II <-> ekliptyczne
# ----------------------------------------------------------------------------
def ecl2eqII( epsilon, lam, beta ):
   """
     współrzędne ekliptyczne (lam,beta)
     -> współrzędne równikowe II (alpha,delta)
     alpha, delta = ecl2eqII( epsilon, lam, beta )
     epsilon - nachylenie równika do ekliptyki
     wszystkie kąty muszą być/są podane w radianach
     na podstawie wzoru cosinusów dla boków trójkąta eliptycznego
     przykład:
     >>> alpha,delta = ecl2eqII( 23.5*d2r, 30.0*d2r, 60.0*d2r )
     >>> round( alpha, 9 )
     -0.261877819
     >>> round( delta, 9 )
     1.10593689
   """
   sd  = np.sin(beta)*np.cos(epsilon) + \
         np.cos(beta)*np.sin(epsilon)*np.sin(lam)
   if ( np.abs((np.abs(sd)-1.0))<eps ):
      print("# <ecl2eqII>: deklinacja zbyt bliska +/-90!" )
      return( 0.0, 0.0 )
   delta  = np.arcsin( sd );
   sa  = (-np.sin(beta)+np.sin(delta)*np.cos(epsilon))\
         /(np.cos(delta)*np.sin(epsilon));
   ca  = np.cos(lam)*np.cos(beta)/np.cos(delta)
   alpha = np.arctan2( sa, ca )
   return( alpha, delta );

def eqII2ecl( epsilon, alpha, delta ):
   """
     współrzędne równikowe II (alpha,delta)
     -> współrzędne ekliptyczne (lam,beta)
     lam, beta = eqII2ecl( epsilon, alpha, delta )
     kąty muszą być podane w radianach
     epsilon - nachylenie równika do ekliptyki
     na podstawie wzoru cosinusów dla boków trójkąta eliptycznego
     przykład:
     >>> lam, beta = eqII2ecl( 23.5*d2r, 30.0*d2r, 60.0*d2r )
     >>> round( lam, 9 ), round( beta, 9 )
     (0.924994875, 0.767738707)
   """
   sb  = np.sin(delta)*np.cos(epsilon)\
        -np.cos(delta)*np.sin(epsilon)*np.sin(alpha)
   if ( np.abs((np.abs(sb)-1.0))<eps ):
      print("# <eqII2ecl>: szerokość ekliptyczna zbyt bliska +/-90!" )
      return( 0.0, 0.0 )
   beta  = np.arcsin( sb );
   sl  = (np.sin(delta)-np.sin(beta)*np.cos(epsilon))\
         /(np.cos(beta)*np.sin(epsilon));
   cl  = np.cos(alpha)*np.cos(delta)/np.cos(beta)
   lam = np.arctan2( sl, cl )
   return( lam, beta );

# ----------------------------------------------------------------------------
# Dzien Juliański, ciągła data Juliańska, średni czas gwiazdowy Greenwich
# ----------------------------------------------------------------------------
def JDN( year, month, day ):
   """
     Dzień Juliański _w południe uniwersalne 12h UTC_
     jdn = JDN( year, month, day )
     Funkcja nie zwraca poprawnej wartości JDN dla dat przed 4 X 1582.
     Dzień Juliański (JDN) to liczba całkowita przypisana do dnia słonecznego
     od południa czasu uniwersalnego Greenwich.
     Dzień juliański JDN=0 to 1. I 4713BC w kalendarzu proleptycznym juliańskim
     (24 XI 4714BC, w kalendarzu proleptycznym gregoriańskim).
     Przejście z kalendarza juliańskiego na gregoriański nastąpiło w X 1582,
     zgodnie z ustaleniami papieża Grzegorza XIII. Dla dat wcześniejszych od
     4. X 1582 lub używany jest kalendarz juliański; dla dat 15. X 1582
     lub później używany jest kalendarz gregoriański; wprowadza to
     10-dniową lukę w datach kalendarzowych, ale jest ciągłość dni juliańskich
     lub dni tygodnia: 4. X 1582 (juliański) to czwartek (JD 2299159.5)
     a 15. X 1582 (gregoriański) to piątek (JD 2299160.5).
     Pominięcie 10 dni w datach kalendarzowych było konieczne ze względu
     na błąd kalendarza juliańskiego narosły przez wiele wieków, wynikający ze
     zbyt częstego wprowadzania lat przestępnych.
     (https://aa.usno.navy.mil/data/JulianDate)
     przykład:
     >>> JDN( 2019, 3, 15 )
     2458558.0
     >>> JDN( 2000, 1, 1 )
     2451545.0
   """
   a = np.floor((14-month)/12.)
   y = year+4800-a
   m = month+12*a-3
   djn = day + np.floor((153*m+2)/5.) + y*365 + \
         np.floor(y/4.) - np.floor(y/100.) + np.floor(y/400.) -32045.0
   return( djn )

def JD( year, month, day, hour, minute, second ):
   """
     Data Juliańska, korzysta z algorytmu numeru Dnia Juliańskiego
     jd = JD( year, month, day, hour, minute, second )
     uwaga:
     funkcja nie zwraca poprawnej wartości JD dla dat przed 4 X 1582
     
     przykład:
     >>> JD( 2019, 3, 13, 0, 0, 0 )
     2458555.5
     >>> JD( 2000, 1,  1, 12,0, 0 )
     2451545.0
   """
   jd = JDN( year, month, day ) + \
        (hour-12.0)/24.0 + minute/1440.0 + second/86400.0
   return( jd )

def GMST( year, month, day, hour, minute, second ):
   """
     średni czas gwiazdowy Greenwich (GMST) dla daty UT1 (UTC) w radianach
     dokładność rzędu 1 sekundy w latach 1900-2100
     GMST( year, month, day, hour, minute, second )
     przykład: rad2hr(GMST(2010,6,10,17,30,30)) -> 10h 46m 20.46s"
     W/g Rocznika Astronomicznego GUGiK         -> 10h 46m 20.39s"
     odchyłka od czasu prawdziwego wynikająca z nutacji jest rzędu 1.1 s
     przykład:
     >>> round(GMST(2010,6,10,17,30,30),12)
     2.820194563323
     >>> printhr(GMST(2010,6,10,17,30,30))
     10h 46m  20.46s
     >>> round( rad2hr(GMST(2010,6,10,17,30,30))[2], 3)
     20.459
   """
   JulianDate  = JD( year, month, day, hour, minute, second )
   D           = JulianDate - JD2000
   time        = np.mod( 18.697374558 + 24.06570982441908*D, 24 )
   return( time*15*np.pi/180.0 )

def GST( year, month, day, hour, minute, second ):
   """
     prawdziwy czas gwiazdowy Greenwich (GAST) dla daty UT1 (UTC) w radianach
     https://aa.usno.navy.mil/faq/docs/GAST.php
     GST( year, month, day, hour, minute, second )
     przykład:
     >>> round( GST( 2019, 3, 13, 0, 0, 0 ), 12 )
     2.972951470668
     >>> printhr(GST( 2019, 3, 13, 0, 0, 0 ))
       11h 21m  21.02s      
     >>> round( rad2hr(GST(2019,3,13,0,0,0))[2], 2)
     21.02
   """
   # argument czasu
   D       = JD( year, month, day, hour, minute, second ) -JD2000
   time    = np.mod( 18.697374558 +24.06570982441908*D, 24 )
   # nachylenie równika do ekliptyki
   epsilon = np.pi*( 23.4393 -0.0000004*D )/180.0
   # węzeł orbity Księżyca
   Omega   = np.pi*( 125.04 -0.052954*D )/180.0
   # średnia długość Słońca
   L       = np.pi*( 280.47 +0.98565*D )/180.0
   # nutacja w długości
   dPsi    = -0.000319*np.sin(Omega) - 0.000024*np.sin(2*L);
   # równanie czasu
   eqeq    = dPsi*np.cos(epsilon)
   #
   time    = time + eqeq
   #
   return( time*15*np.pi/180.0 )

def JDtoDate(jd):
    """
      zamienia datę juliańską na datę kalendarzową
      algoritm 'Practical Astronomy with your Calculator or Spreadsheet',
        4th ed., Duffet-Smith and Zwart, 2011.
     JDtoDate(jd)
     parametery:
     jd : float
        Julian Day
     zwraca
     year : int
        Year as integer. Years preceding 1 A.D. should be 0 or negative.
        The year before 1 A.D. is 0, 10 B.C. is year -9.
     month : int
        Month as integer, Jan = 1, Feb. = 2, etc.
     day : float
        Day, may contain fractional part.
     
     przykład:
     ---------
     convert Julian Day 2446113.75 to year, month, and day.
     >>> JDtoDate(2446113.75)
     (1985.0, 2.0, 17.25)
     >>> JDtoDate(2451545.0)
     (2000.0, 1.0, 1.5)
     >>> JD( 2000.0, 1.0, 1, 12, 0, 0 )
     2451545.0
     >>> asa.JDtoDate( 0 )
     (-4712.0, 1.0, 1.5)
    """
    jd = jd + 0.5
    F, I = np.modf(jd)
    I = int(I)
    A = np.trunc((I - 1867216.25)/36524.25)
    if I > 2299160:
        B = I + 1 + A - np.trunc(A / 4.)
    else:
        B = I
    C = B + 1524
    D = np.trunc((C - 122.1) / 365.25)
    E = np.trunc(365.25 * D)
    G = np.trunc((C - E) / 30.6001)
    day = C - E + F - np.trunc(30.6001 * G)
    if G < 13.5:
        month = int(G - 1)
    else:
        month = int(G - 13)
    if month > 2.5:
        year = int(D - 4716)
    else:
        year = int(D - 4715)
    return year, month, day

# ----------------------------------------------------------------------------
# Testy modułu
# ----------------------------------------------------------------------------
def opis(fun):
  """
    skrócona forma wywołania opisu funkcji w module
  """
  print("\nOpis funkcji", fun )
  print(inspect.getdoc(fun))

def doc(fun):
  """
    skrócona forma wywołania dokumentacji funkcji w module
  """
  print( inspect.getdoc(fun) )

# ----------------------------------------------------------------------------
# Konwersja elementów wektorowych
# ----------------------------------------------------------------------------
def sph2cart( r, lon, phi ):
   """
     konwersja współrzędnych biegunowych [r,lon,phi] na kartezjańske [x,y,z]
     x,y,z = sph2cart( r, lon, phi )
     wersja skalarna
     przykład:
     >>> phi=deg2rad(35,15,51.80)
     >>> lon=deg2rad(45,0,0)
     >>> r=1.7320508075688772
     >>> sph2cart( r, lon, phi )
     (1.0000000097973685, 1.0000000097973683, 0.9999999804052628)
   """
   return( r*np.cos(phi)*np.cos(lon), r*np.cos(phi)*np.sin(lon), r*np.sin(phi))

def cart2sph( x, y, z ):
   """
     konwersja współrzędnych kartezjańskich [x,y,z] na biegunowe [r,lambda,phi]
     r, lam, phi = cart2sph( x, y, z )
     wersja skalarna
     przykład:
     >>> r, lon, phi = cart2sph( 1,1,1 )
     >>> print(r)
     1.7320508075688772
     >>> printdeg(lon)
     45  0' 0.00"
     >>> printdeg(phi)
     35 15' 51.80"
   """
   rho  = np.sqrt( x*x + y*y + z*z )
   lon  = np.mod( np.arctan2( y/rho, x/rho ), 2*np.pi )
   phi  = np.mod( np.arcsin( z/rho ), 2*np.pi )
   return ( rho, lon, phi )

def sph2cartv( r, lon, phi ):
   """
     konwersja współrzędnych biegunowych [r,lon,phi] na kartezjańskie [x,y,z]
     kąt phi ma sens szerokości, lon - długości, kąty podane w radianach
     v = sph2cartv( r, lon, phi )
     wersja wektorowa
   """
   s = r*np.array( [\
       np.cos(phi)*np.cos(lon), np.cos(phi)*np.sin(lon), np.sin(phi) \
       ] )
   return( s )

def cart2sphv( r ):
   """
     konwersja współrzędnych kartezjańskich [x,y,z] na biegunowe [r,lambda,phi]
     wersja wektorowa x = r[1], y = r[2], z = r[3]
     v = cart2sphv( r )
     przykład:
     >>> phi=deg2rad(35,15,51.80)
     >>> lon=deg2rad(45,0,0)
     >>> r=1.7320508075688772
     >>> sph2cartv( r, lon, phi )
     array([1.00000001, 1.00000001, 0.99999998])
   """
   rho  = np.sqrt( r[0]*r[0] + r[1]*r[1] + r[2]*r[2] )
   lon  = np.mod( np.arctan2( r[1]/rho, r[0]/rho ), 2*np.pi )
   phi  = np.mod( np.arcsin( r[2]/rho ), 2*np.pi )
   return( rho, lon, phi )

# ----------------------------------------------------------------------------
# macierze obrótow elementarnych
# ----------------------------------------------------------------------------
def A1( phi ):
   """
     macierz obrotu elementarnego wokół osi 1 o kąt dodatni phi
     dla obrotów o kat phi wokół osi 1 -> x
     przykład:
     >>> A1( deg2rad(30,0,0) )
     array([[ 1.       ,  0.       ,  0.       ],
              [ 0.       ,  0.8660254,  0.5      ],
              [ 0.       , -0.5      ,  0.8660254]])
   """
   A = np.array( [ [     1,       0,           0        ],\
                   [     0, +np.cos(phi), +np.sin(phi)  ],\
                   [     0, -np.sin(phi), +np.cos(phi)  ] ] )
   return( A )

def A2( theta ):
   """
     macierz obrotu elementarnego wokół osi 2 o kat dodatni phi
     dla obrotow o kat phi wokół osi 2 -> y
     przykład:
     >>> A2( deg2rad(30,0,0) )
     array([[ 0.8660254,  0.       , -0.5      ],
            [ 0.       ,  1.       ,  0.       ],
            [ 0.5      ,  0.       ,  0.8660254]])
   """
   A = np.array([ [ +np.cos(theta), 0, -np.sin(theta) ],\
                  [     0,          1,        0       ],\
                  [ +np.sin(theta), 0, +np.cos(theta) ] ] )
   return( A )

def A3( psi ):
   """
     macierz obrotu elementarnego wokół osi 3 o kąt dodatni phi
     dla obrotów o kat phi wokół osi 3 -> z
     przykład:
     >>> A3( deg2rad(30,0,0) )
     array([[ 0.8660254,  0.5      ,  0.       ],
            [-0.5      ,  0.8660254,  0.       ],
            [ 0.       ,  0.       ,  1.       ]])
   """
   A = np.array( [ [ +np.cos(psi), +np.sin(psi), 0 ],\
                   [ -np.sin(psi), +np.cos(psi), 0 ],\
                   [       0,           0,       1 ] ] )
   return( A )

# ----------------------------------------------------------------------------
# Macierz precesji ogólnej
# ----------------------------------------------------------------------------
def precess( dt ):
   """
     oblicza macierz precesji ogólnej P w/g teorii IAU 1986
                        r( t0+dt ) = P * r(t0)
     czas dt jest liczony w DNIACH juliańskich względem epoki J2000.0
     precess( dt )
     przykład:
     >>> precess( 22*365.25 )
     array([[ 9.99985612e-01, -4.91988179e-03, -2.13775978e-03],
            [ 4.91988179e-03,  9.99987897e-01, -5.25899848e-06],
            [ 2.13775979e-03, -5.25860262e-06,  9.99997715e-01]])
   """
   t      = dt/cy
   # kąty precesji ksiA, zA, thetaA
   ksiA   = ( 0.6406161*t + 0.0000839*t*t + 0.0000050*t*t*t )*d2r
   zA     = ( 0.6406161*t + 0.0003041*t*t + 0.0000005*t*t*t )*d2r
   thetaA = ( 0.5567530*t - 0.0000118*t*t - 0.0000011*t*t*t )*d2r
   # macierz precesji ogolnej
   P      = np.dot( A3(-zA), np.dot( A2(thetaA), A3(-ksiA) ) )
   return( P )

# ----------------------------------------------------------------------------
# Długość ekliptyczna Słońca i równanie czasu
# ----------------------------------------------------------------------------
def SunLongitude( year, month, day, hour, min, sec ):
   """
     długość ekliptyczna Słońca na podstawie Almagest, R. Fitzpatrick
     szczegóły takiże w wykładzie IX
     dokładność 0.1'-0.7' w latach 1800:2050
     uwzględniona precesja punktu Barana, precesja peryhelium  i aberracja
     testy stałych teorii
     ruch średni w długości po uwzględnieniu precesji:
     >>> (360.0)/365.242189
     0.9856473617838273
     wynik ten sam, gdy do ruchu średniego Ziemi (JPL) dodamy tempo precesji:
     >>> 35999.37244981/36525+3.8246e-5
     0.9856473479797399     
     lub w wersji z wyznaczeniem tempa precesji     
     >>> 35999.37244981/36525.0-0.32327364/36525.0
     0.9856002512298425
     poprawka aberracyjna do długości średniej Ziemi w epoce początkowej
     >>> 100.46457166-20.5/3600.0
     100.45887721555556
   """         
   # stałe modelu orbitalnego Słońca w/g danych JPL
   deg2rad = np.pi/180.0
   t0      = 2451545.0    # epoka początkowa [JD]
   ecc     = 0.01671123   # mimośród obity Ziemi w epoce t0  
   omega0  = 102.93768193 # długość preicentrum [deg] w epoce t0   
   kappa   = 20.5/3600.0  # poprawka aberracyjna [deg]
   psi     = 3.8246e-5    # roczne tempo precesji [deg]
                          # dw/dt precesja peryhelium [deg/dzień]
   nomega  = 0.32327364/36525.0
                          # dL/dt ruch średni w długości                        
   n       = 35999.37244981/36525.0
   nL      = n + psi      # ruch średni w długości z precesją [deg/dzień]                 
   nM      = n - nomega   # ruch średni z precesją peryhelium [deg/dzień]
                          # długość średnia epoki + popr. aberracji [deg]
   lambda0 = 100.46457166-kappa
                          # anomalia średnia w epoce t0
   M0      = lambda0 - omega0
   #
   # Algorytm obliczenia długości ekliptycznej Słońca
   #
   # data juliańska momentu obserwacji
   t       = JD( year, month, day, hour, min, sec )   
   # długość średnia Ziemi w epoce z aberracją, ruch średni uwzględnia precesję
   meanL   = np.mod( lambda0 + nL*(t-t0), 360.0 )
   # anomalia śednia w epoce obserwacji, uwzględnia precesję peryhelium
   M       = np.mod( M0 + nM*(t-t0), 360.0 )*deg2rad
   # poprawka ze względu na równanie środka, wyznaczające anomalię prawdziwą
   q       = 2.0*ecc*np.sin(M) + (5.0/4.0)*ecc*ecc*np.sin(2.0*M)
   #
   # dokładniejsza wersja równania środka, tutaj bez wpływu na końcowy wynik
   #
   #q      = (2.0*ecc-ecc**3/4.0)*np.sin(M) + \
   #         (5.0/4.0)*ecc*ecc*np.sin(2.0*M) + \
   #         (13.0/12.0)*ecc*ecc*ecc*np.sin(3.0*M)
   #
   # długość ekliptyczna Słońca w układzie równika i ekliptyki daty
   return( np.mod( meanL + q/deg2rad +180.0, 360.0 )*deg2rad );

def SunLambda( year, month, day, hour, min, sec ):
   """
     długość ekliptyczna Słońca [rad], dokładność 1' w 1800:2050
     uwzględniona precesja, precesja peryhelium i aberracja
     algorytm z Astronomical Almanac, sekcja C24  
     przykład:
     >>> lsun1 = SunLongitude( 2024, 6, 4, 10, 22, 22 )
     >>> lsun2 = SunLambda( 2024, 6, 4, 10, 22, 22 )
     >>> print180( lsun1-lsun2 )
     -0  0' 0.08"
   """         
   # interwał od daty juliańskiej t0 = 2451545.0
   dt    = JD( year, month, day, hour, min, sec ) - 2451545.0
   # długość średnia Słońca, poprawiona na precesję i aberrację
   meanL = np.mod( 280.460 + 0.9856474*dt, 360.0 )
   # anomalia śednia w epoce obserwacji z precesją peryhelium
   M     = np.mod( 357.528 + 0.9856003*dt, 360.0 )*np.pi/180.0
   # poprawka na równanie środka wyznaczające anomalię prawdziwą
   q     = 1.915*np.sin(M) + 0.020*np.sin(2.0*M)
   # długość ekliptyczna Słońca w układzie równika-ekliptyki daty
   return( np.mod( meanL + q, 360.0 )*np.pi/180.0 );

def SunLongitudeJD( jd ):
   """
     długość ekliptyczna Słońca z argumentem Daty Juliańskiej
     dokładność 0.1'-0.7' w latach 1800:2050
     uwzględniona precesja punktu Barana, precesja peryhelium obity i aberracja
     algorytm z Almagest, R. Fitzpatrick
     ruch średni w długości po uwzględnieniu precesji:
     >>> (360.0)/365.242189
     0.9856473617838273
     wynik ten sam, gdy do ruchu średniego Ziemi (JPL) dodamy tempo precesji:
     >>> 35999.37244981/36525+3.8246e-5
     0.9856473479797399     
     lub w wersji z wyznaczeniem tempa precesji     
     >>> 35999.37244981/36525.0-0.32327364/36525.0
     0.9856002512298425
     poprawka aberracyjna do długości średniej Ziemi w epoce początkowej
     >>> 100.46457166-20.5/3600.0
     100.45887721555556
   """         
   # stałe modelu orbitalnego Słońca w/g danych JPL
   deg2rad = np.pi/180.0
   t0      = 2451545.0    # epoka początkowa [JD]
   ecc     = 0.01671123   # mimośród obity Ziemi w epoce t0  
   omega0  = 102.93768193 # długość pericentrum [deg] w epoce t0   
   kappa   = 20.5/3600.0  # poprawka aberracyjna [deg]
   psi     = 3.8246e-5    # roczne tempo precesji lunisolarnej [deg]
                          # dw/dt precesja peryhelium [deg/dzień]
   nomega  = 0.32327364/36525.0
                          # dL/dt ruch średni w długości                        
   n       = 35999.37244981/36525.0
   nL      = n + psi      # ruch średni w długości z precesją [deg/dzień]                 
   nM      = n - nomega   # ruch średni z precesją peryhelium [deg/dzień]
                          # długość średnia epoki + popr. aberracji [deg]
   lambda0 = 100.46457166-kappa
                          # anomalia średnia w epoce t0
   M0      = lambda0 - omega0
   #
   # Algorytm obliczenia długości ekliptycznej Słońca
   #
   # data juliańska momentu obserwacji
   t       = jd
   # długość średnia Ziemi w epoce z aberracją, ruch średni uwzględnia precesję
   meanL   = np.mod( lambda0 + nL*(t-t0), 360.0 )
   # anomalia śednia w epoce obserwacji, uwzględnia precesję peryhelium
   M       = np.mod( M0 + nM*(t-t0), 360.0 )*deg2rad
   # poprawka ze względu na równanie środka, wyznaczające anomalię prawdziwą
   q       = 2.0*ecc*np.sin(M) + (5.0/4.0)*ecc*ecc*np.sin(2.0*M)
   #
   # dokładniejsza wersja równania środka, tutaj bez wpływu na końcowy wynik
   #
   #q      = (2.0*ecc-ecc**3/4.0)*np.sin(M) + \
   #         (5.0/4.0)*ecc*ecc*np.sin(2.0*M) + \
   #         (13.0/12.0)*ecc*ecc*ecc*np.sin(3.0*M)
   #
   # długość ekliptyczna Słońca w układzie równika i ekliptyki daty
   return( np.mod( meanL + q/deg2rad +180.0, 360.0 )*deg2rad );

def DeltaE( year, month, day, hour, min, sec ):
   """
     równanie czasu w oparciu o długość ekliptyczną Słońca [min]
              ΔE = T_true - T_mean = α_mean - α_true
     długość ekliptyczna Słońca obliczana w/g Almagest lub AA
     przykład:
     >>> lsun = SunLongitude( 2024, 6, 4, 10, 22, 22 )
     >>> printdeg( lsun, "λ☉ w dniu 6.4.2024, godz. 10:22:22 UTC)")
     λ☉ w dniu 6.4.2024, godz. 10:22:22 UTC) =  74  0' 17.54"
     >>> jt = JD( 2024, 6, 4, 10, 22, 22 )
     2460465.932199074
     >>> lsun = SunLongitudeJD( jt )
     >>> printdeg( lsun, "λ☉ w dniu 6.4.2024, godz. 10:22:22 UTC)")
     λ☉ w dniu 6.4.2024, godz. 10:22:22 UTC) =  74  0' 17.63"
   """    
   # stałe modelu orbitalnego Słońca w/g danych JPL
   ecc     = 0.01671123   # mimośród obity Ziemi   
   lambda0 = 280.458      # długość średnia epoki z poprawką aberracyjną [deg]
   M0      = 357.588      # anomalia średnia epoki   
   eps     = deg2rad( 23, 26, 21.45 )
   # współrzędne kątowe Słońca: długość średnia
   lsun = SunLongitude( year, month, day, hour, min, sec );
   # i anomalia średnia, M = lsun + M0-lambda0, dla wyliczenia q
   M    = np.mod( lsun + (M0-lambda0)*np.pi/180, 2*np.pi );
   # rektascensja RA Słońca ,,prawdziwego''
   asun = np.arctan( np.cos(eps)*np.tan(lsun) )
   if ( lsun >= 3*np.pi/2.0 ):
     alphasun = 2*np.pi+asun
   if ( lsun <= np.pi/2.0 ):
     alphasun = asun
   if ( (lsun > np.pi/2.0) & ( lsun < 3*np.pi/2.0) ):   
     alphasun = asun + np.pi
   # równanie czasu: mean alpha - alpha = lambda - arctan( alpha ) + q
   dalpha = lsun -alphasun -2*ecc*np.sin(M)
   return( (dalpha*180.0/np.pi)/15.0*60 );   

def refract( zobs, Hobs=0, P=1013, T=0):
   """
     oblicza refrakcję atmosferyczną w kącie odległości zenitalnej
     dla długości fal optycznych, wynik w radianach
     dodatek wyjaśniający do Astronomical Almanach, str. 144.
     Bennett, G. (1982), Journal of Navigation (Royal Institute) 35,
     s. 255-259.
     zobs: odległość zenitalna z observacji w radianach
     Hobs: wysokość nad poziomem morza w m
     na podstawie NOVAS, V1.0/06-98/JAB (USNO/AA)
     funkcja może być używana do planowania obserwacji lub
     sterowania teleskopem, ale nie powinna być stosowana do redukcji
     precyzyjnych obserwacji.
     przykład:
     >>> asa.printdeg( asa.refract( 10*np.pi/180., 0 ) )
     0  0' 10.90"
     >>> asa.printdeg( asa.refract( 40*np.pi/180., 0 ) )
     0  0' 52.06"
     >>> asa.printdeg( asa.refract( 80*np.pi/180., 0 ) )
     0  5' 36.10"
     >>> asa.printdeg( asa.refract( 90*np.pi/180., 0 ) )
     0 35' 49.32"   
   """
   # 's' to przybliżona wysokość atmosfery w metrach
   s = 9.1e3
   zo = zobs * 180/np.pi
   # refrakcja tylko dla odległości zenitalnych pomiędzy 0.1 a 91 stopni
   if ((zobs < 0.1) | (zobs > 91.0)):
       refr = 0.0
   # jeśli dostępne są dane pogodowe z obserwacji, to użyj ich,
   # w przeciwnym razie użyj szacunków warunków średnich
   p = P
   t = T
   if ( (P==1010) & (T==0) ):
       p = P * np.exp( Hobs / s )
   h = 90.0 - zo
   r = 0.016667 / np.tan( (h + 7.31/(h+4.4) ) * np.pi/180.)
   refr = r * (0.28 * p / (t + 273.0)) * np.pi/180.
   return(  refr )

# ----------------------------------------------------------------------------
# Długość ekliptyczna Księżyca i fazy Księżyca
# ----------------------------------------------------------------------------
def MoonLongitude( jd ):
  """
    długość ekliptyczna Księżyca dla daty juliańskiej JD przedstawiony
    jako średnia orbita + perturbacje harmoniczne
    algorytm z książki Jana Meeusa
    Chapront ELP-2000/82
    pominięty wpływ Wenus, Jowisza i spłaszczenia Ziemi
    dokładność 0.1 stopnia (6') w interwale +/-100 lat od J2000.0
  """
  T = (jd - 2451545.0)/36525.0
  # elementy orbity Księżyca (Meeus) względem J2000
  # średnia długość ekliptyczna  
  Lp = np.deg2rad( np.mod(
        218.3164477 + 481267.88123421*T
        - 0.0015786*T**2 + T**3/538841 - T**4/65194000, 360) )
  # średnia elongacja (różnica średnich długości Słońa i Księżyca)
  D = np.deg2rad( np.mod( 297.8501921 + 445267.1114034*T
        - 0.0018819*T**2 + T**3/545868 - T**4/113065000, 360) )
  # anomalia średnia Słońca
  M = np.deg2rad( np.mod( 357.5291092 + 35999.0502909*T
        - 0.0001536*T**2 + T**3/24490000, 360) )
  # anomalia średnia Księżyca
  Mp = np.deg2rad( np.mod( 134.9633964 + 477198.8675055*T
        + 0.0087414*T**2 + T**3/69699 - T**4/14712000, 360) )
  # argument szerokości (separacja od węzła wstępującego)
  F = np.deg2rad( np.mod( 93.2720950 + 483202.0175233*T
        - 0.0036539*T**2 - T**3/3526000 + T**4/863310000, 360) )
  # ---- perturbacje w długości ekliptycznej (deg)
  dlam = ( 6.289*np.sin( Mp )     # równanie środka
        + 1.274*np.sin( 2*D - Mp )  # ewekcja, słynna w XVII w.
        + 0.658*np.sin( 2*D )       # wariacja
        + 0.214*np.sin( 2*Mp )      # + mniejsze efekty
        + 0.110*np.sin( D )
  )
   
  lam = Lp + np.deg2rad( dlam )
   
  # ---- perturbacje w szerokości ekliptycznej (deg)
  beta = np.deg2rad( 5.128*np.sin( F ) # od inklinacji 5.128 deg
        + 0.280*np.sin( Mp + F )
        + 0.277*np.sin( Mp - F )
        + 0.173*np.sin( 2*D - F )
        + 0.055*np.sin( 2*D + F - Mp )
        + 0.046*np.sin( 2*D - F - Mp )
        + 0.033*np.sin( 2*D + F )
        + 0.017*np.sin( 2*Mp + F )
  )
  lam = np.mod( lam, 2*np.pi )
  return( lam, beta )

def NewMoonSimple( jd0 ):
  """
    moment nowiu najbliższy zadanej dacie juliańskiej
    różnica ruchów średnich Księżyca i Słońca
    dλ_M/dt − dλ⊙/dt ≈ 12.19∘/day
    przy czym
    d λ_M/dt ≈  0.986∘/d
    dλ⊙/dt ≈ 12.190∘/d
    ze znanej różnicy długości Δλ=λ_M−λ⊙ mamy
    Δt=−Δλ/12.190∘/d
  """
  jd = jd0
  for i in range( 5 ):
        lam_moon, _ = MoonLongitude( jd )
        lam_sun = SunLongitudeJD( jd )
        dlam = wrap( lam_moon - lam_sun )
        # poprawka iteracyjna vs Δ ruchów średnich 12.19 deg/day
        jd -= dlam / np.deg2rad( 12.19075 )
  return( jd )

def NewMoon( jd0 ):
  """
    moment nowiu w oparciu o warunek koniunkcji ekliptycznej
    λ_Moon(t) = λ_Sun(t) czyli f(t) = λ_Moon(t) - λ_Sun(t) = 0
    algorytm oparty o metodę Newtona-Raphsona
  """
  jd = jd0
  for i in range( 6 ):
     lam_moon, _ = MoonLongitude( jd )
     lam_sun = SunLongitudeJD( jd )
     f = wrap( np.deg2rad(lam_moon - lam_sun) )
     # pochodna numeryczna i algorytm Newtona-Raphsona
     dt = 0.001
     lam_moon2,_ = MoonLongitude( jd + dt )
     lam_sun2 = SunLongitudeJD( jd + dt )
     f2 = wrap( np.deg2rad( lam_moon2 - lam_sun2) )
     df = (f2 - f)/dt
     jd -= f/df
  return( jd )    

def MoonPhaseTime( jd0, phase ):
  """
    fazy Księżyca dla argumentu phase w stopniach:
      nów phase => 0 deg
      pierwsza i ostatnia kwadra => 90 deg i 270 deg
      pełnia => 180 deg
    podlegają okresowi synodycznemu
    Psyn​ = 360∘/12.19∘​ ≈ 29.53d  
    algorytm Jana Meeusa
  """
  jd = jd0
  phi = np.deg2rad( phase )
  for i in range(5):
    lam_moon,_ = MoonLongitude( jd )
    lam_sun = SunLongitudeJD( jd )
    f = wrap( lam_moon - lam_sun - phi )
    # średnia prędkość (zmiana) fazy
    jd -= f / np.deg2rad( 12.19075 )
  return( jd )   

def test():
  """
    funkcja testowa modułu astrobox
  """
  print("\n" + "="*70)
  print("TESTOWANIE MODUŁU ASTROBOX")
  print("="*70)
  
  print("\n--- Test 1: Lista funkcji modułu ---")
  print("Lista funkcji: dir(asa) ")
  print("\nPrzykład opisu funkcji JDN:")
  print(inspect.getdoc(JDN))
  
  print("\n--- Test 2: Konwersje kątów i daty juliańskie ---")
  DEC = deg2rad( 89, 15, 40.8 )
  printdeg( DEC, 'deklinacja Polaris' )
  print("deklinacja Polaris, bez nazwy:")
  printdeg( DEC )
  
  print("\nData juliańska:")
  jd_test = JD( 2020, 3, 18, 14, 10, 10 )
  print(f"JD(2020, 3, 18, 14, 10, 10) = {jd_test}")
  
  print("\n--- Test 3: Czas gwiazdowy ---")
  gmst = GMST( 2019, 3, 18, 14, 0, 0 )
  printhr( gmst, 'średni czas gwiazdowy Greenwich   ' )
  gst = GST( 2019, 3, 18, 14, 0, 0 )
  printhr( gst,  'prawdziwy czas gwiazdowy Greenwich' )
  
  print("\n--- Test 4: Przemiana ekliptyczne -> równikowe II ---")
  alpha,delta = ecl2eqII( 23.5*d2r, 30.0*d2r, 23.5*d2r )
  print( f"alpha = {alpha:.6f}, delta = {delta:.6f}" )
  
  print("\n--- Test 5: Długość ekliptyczna Słońca w równonocy (model Almagest) ---")
  print("Daty z rocznika astronomicznego:")
  printdeg( SunLongitude( 1896,3,20,  2,46,0 ), "equinox 1896 02:46, λ☉ = " );
  printdeg( SunLongitude( 2010,3,20, 17,32,0 ), "equinox 2010 17:32, λ☉ = " );
  printdeg( SunLongitude( 2015,3,20, 22,45,0 ), "equinox 2015 22:45, λ☉ = " );
  printdeg( SunLongitude( 2016,3,20,  4,30,0 ), "equinox 2016 04:30, λ☉ = " );
  printdeg( SunLongitude( 2019,3,20, 21,58,0 ), "equinox 2019 21:58, λ☉ = " );
  printdeg( SunLongitude( 2020,3,20,  3,49,0 ), "equinox 2020 03:49, λ☉ = " );
  
  print("\n--- Test 6: Długość ekliptyczna Słońca (model AA) ---")
  print("Daty z rocznika, Astronomical Almanac:")
  printdeg( SunLambda( 1896, 3, 20,  2, 46,41 ), "equinox 1896 02:46, λ☉ = " );
  printdeg( SunLambda( 2010, 3, 20, 17, 32, 0 ), "equinox 2010 17:32, λ☉ = " );
  printdeg( SunLambda( 2015, 3, 20, 22, 45, 0 ), "equinox 2015 22:45, λ☉ = " );
  printdeg( SunLambda( 2016, 3, 20,  4, 30, 0 ), "equinox 2016 04:30, λ☉ = " );
  printdeg( SunLambda( 2019, 3, 20, 21, 58, 0 ), "equinox 2019 21:58, λ☉ = " );
  printdeg( SunLambda( 2020, 3, 20,  3, 49, 0 ), "equinox 2020 03:49, λ☉ = " );  
  
  print("\n--- Test 7: Słońce dla momentów z AA, układ średni daty ---")
  printdeg( SunLongitude( 2000,3,20,  7,35,0 ), "equinox 2000 20.III  7:35 ");
  printdeg( SunLongitude( 2000,9,22, 17,28,0 ), "equinox 2000 22.IX  17:28 ");
  printdeg( SunLongitude( 2010,3,20, 17,32,0 ), "equinox 2010 20.III 17:32 ");
  printdeg( SunLongitude( 2010,9,23,   3,9,0 ), "equinox 2010 23.IX   3:09 ");
  printdeg( SunLongitude( 2015,3,20, 22,45,0 ), "equinox 2015 20.III 22:45 ");
  printdeg( SunLongitude( 2020,3,20,  3,49,0 ), "equinox 2020 20.III 03:49 ");   
  printdeg( SunLongitude( 2020,9,22, 13,31,0 ), "equinox 2020 22.IX  13:31 ");
  printdeg( SunLongitude( 2024,3,20,   3,6,0 ), "equinox 2024 20.III 03:06 ");   
  printdeg( SunLongitude( 2024,9,22, 12,44,0 ), "equinox 2024 22.IX  12:44)");
  
  print("\n--- Test 8: Równanie czasu w 2024 roku (co tydzień) ---")
  print("Odchyłka 1' kątowej = 1/30 kątowych rozmiarów tarczy Słońca")
  print("Równania czasu w kolejnych tygodniach 2024, podane w minutach:")
  year = 2024
  hr   = 12
  mday = ( 31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31 )
  total_days = sum(mday)
  date = np.zeros( total_days )
  de   = np.zeros( total_days )
  k    = 0
  print(f"{'Dzień roku':<12} {'Równanie czasu':<18} {'Data':<15}")
  print("-" * 45)
  for i in range( 0, 12 ):
    for j in range( 0, mday[i], 7 ):
       date[k] = JD( year, i+1, j+1, hr, 0, 0 ) - JD( year, 1, 1, hr, 0, 0 )
       de[k]   = DeltaE( year, i+1, j+1, hr, 0, 0 )
       print(f"{date[k]:<12.2f} {DeltaE( year, i+1, j+1, hr, 0, 0 ):<18.3f} {year}-{i+1:02d}-{j+1:02d}")
       k = k+1
       
  print("\n--- Test 9: Konwersje współrzędnych ---")
  print("Sferyczne -> Kartezjańskie:")
  phi_test = deg2rad(35,15,51.80)
  lon_test = deg2rad(45,0,0)
  r_test = 1.7320508075688772
  cart = sph2cart( r_test, lon_test, phi_test )
  print(f"sph2cart({r_test}, {lon_test:.6f}, {phi_test:.6f}) = {cart}")
  
  print("\nKartezjańskie -> Sferyczne:")
  sph = cart2sph( 1, 1, 1 )
  print(f"cart2sph(1, 1, 1) -> r={sph[0]:.6f}, lon={sph[1]:.6f}, phi={sph[2]:.6f}")
  
  print("\n--- Test 10: Macierze obrotów ---")
  A1_test = A1( deg2rad(30,0,0) )
  print(f"A1(30) =\n{A1_test}")
  
  print("\n--- Test 11: Precesja ---")
  P = precess( 22*365.25 )
  print(f"Macierz precesji dla dt=22 lat juliańskich:")
  print(P)
  
  print("\n--- Test 12: Księżyc - długość ekliptyczna ---")
  jd_moon = 2459580.5  # Przykładowa data juliańska
  lam_moon, beta_moon = MoonLongitude( jd_moon )
  print(f"MoonLongitude({jd_moon}) -> λ={np.degrees(lam_moon):.4f}, β={np.degrees(beta_moon):.4f}")
  
  print("\n--- Test 13: Fazy Księżyca ---")
  # Nowy Księżyc - styczeń 2024
  jd_new_moon = NewMoon( JD(2024, 1, 1, 0, 0, 0) )
  date_new = JDtoDate(jd_new_moon)
  print(f"Nowy Księżyc najbliższy 2024-01-01: JD={jd_new_moon:.4f}")
  print(f"Data: {int(date_new[0])}-{int(date_new[1]):02d}-{date_new[2]:.2f}")
  
  print("\n--- Test 14: Refrakcja atmosferyczna ---")
  z1 = 10 * np.pi / 180.
  z2 = 40 * np.pi / 180.
  z3 = 80 * np.pi / 180.
  print(f"refract(10) = {np.degrees(refract(z1, 0))*3600:.2f}\"")
  print(f"refract(40) = {np.degrees(refract(z2, 0))*3600:.2f}\"")
  print(f"refract(80) = {np.degrees(refract(z3, 0))*3600:.2f}\"")
  
  print("\n--- Test 15: Konwersja daty juliańskiej ---")
  jd_test = 2451545.0
  cal_date = JDtoDate(jd_test)
  print(f"JDtoDate({jd_test}) = {cal_date}")
  
  jd_test2 = 2446113.75
  cal_date2 = JDtoDate(jd_test2)
  print(f"JDtoDate({jd_test2}) = {cal_date2}")
  
  print("\n--- Test 16: Konwersja dni na godziny ---")
  dt_test = 1.5
  result = day2hr(dt_test)
  print(f"day2hr({dt_test}) = {result}")
  
  dt_test2 = 0.123456
  result2 = day2hr(dt_test2)
  print(f"day2hr({dt_test2}) = {result2}")
  
  print("\n" + "="*70)
  print("TESTY ZAKOŃCZONE POMYŚLNIE")
  print("="*70 + "\n")

if __name__ == "__main__":
  print("Moduł astrobox, Krzysztof Goździewski, 2009-2024, ver 6.04.2024")
  import doctest
  import inspect
  print("\nTest zgodności wybranych funkcji modułu astrobox.py\n")
  doctest.testmod(verbose=True)
  test()
  print("\n")
# print("Moduł astrobox, Krzysztof Goździewski, 2009-2024, ver 13.04.2024")
# doc(intro)
# ----------------------------------------------------------------------------
# Notatki
# ----------------------------------------------------------------------------
#
# Definicja i algorytm obliczania czasu gwiazdowego (Astronomical Almanac)
#
# https://aa.usno.navy.mil/faq/docs/GAST.php
#
# Sidereal time is a system of timekeeping based on the rotation of the Earth
# with respect to the fixed stars in the sky.  More specifically, it is the
# measure of the hour angle of the vernal equinox.  If the hour angle is
# measured with respect to the true equinox, apparent sidereal time is being
# measured.  If the hour angle is measured with respect to the mean equinox,
# mean sidereal time is being measured.  When the measurements are made with
# respect to the meridian at Greenwich, the times are referred to as
# Greenwich mean sidereal time (GMST) and Greenwich apparent sidereal time
# (GAST).
#
# Given below is a simple algorithm for computing apparent sidereal time to
# an accuracy of about 0.1 second, equivalent to about 1.5 arcseconds on the
# sky.  The input time required by the algorithm is represented as a Julian
# date (Julian dates can be used to determine Universal Time.)
# Let JD be the Julian date of the time of interest.  Let JD0 be the Julian
# date of the previous midnight (0h) UT (the value of JD0 will end in .5
# exactly), and let H be the hours of UT elapsed since that time.  Thus we
# have JD = JD0 + H/24.
# For both of these Julian dates, compute the number of days and fraction (+
# or -) from 2000 January 1, 12h UT, Julian date 2451545.0:
# D = JD - 2451545.0
# D0 = JD0 - 2451545.0
# Then the Greenwich mean sidereal time in hours is
# GMST = 6.697374558 + 0.06570982441908 D0 + 1.00273790935 H + 0.000026 T2
# where T = D/36525 is the number of centuries since the year 2000; thus the
# last term can be omitted in most applications.  It will be necessary to
# reduce GMST to the range 0h to 24h.  Setting H = 0 in the above formula
# yields the Greenwich mean sidereal time at 0h UT, which is tabulated in The
# Astronomical Almanac.
# The following alternative formula can be used with a loss of precision of
# 0.1 second per century:
# GMST = 18.697374558 + 24.06570982441908 D
# where, as above, GMST must be reduced to the range 0h to 24h.  The
# equations for GMST given above are adapted from those given in Appendix A
# of USNO Circular No.  163 (1981).
# The Greenwich apparent sidereal time is obtained by adding a correction to
# the Greenwich mean sidereal time computed above.  The correction term is
# called the nutation in right ascension or the equation of the equinoxes.
# Thus,
# GAST = GMST + eqeq.
# The equation of the equinoxes is given as eqeq = Δψ cos ε where Δψ, the
# nutation in longitude, is given in hours approximately by
# Δψ ≈ -0.000319 sin Ω - 0.000024 sin 2L
# with Ω, the Longitude of the ascending node of the Moon, given as
# Ω = 125.04 - 0.052954 D,
# and L, the Mean Longitude of the Sun, given as
# L = 280.47 + 0.98565 D.
# ε is the obliquity and is given as ε = 23.4393 - 0.0000004 D.
# The above expressions for Ω, L, and ε are all expressed in degrees.
# The mean or apparent sidereal time locally is found by obtaining the local
# longitude in degrees, converting it to hours by dividing by 15, and then
# adding it to or subtracting it from the Greenwich time depending on whether
# the local position is east (add) or west (subtract) of Greenwich.
# If you need apparent sidereal time to better than 0.1 second accuracy on a
# regular basis, consider using the Multiyear Interactive Computer Almanac,
# MICA.  MICA provides very accurate almanac data in tabular form for a range
# of years.
#
# NOTES ON ACCURACY
# The maximum error resulting from the use of the above formulas for sidereal
# time over the period 2000-2100 is 0.432 seconds; the RMS error is 0.01512
# seconds.  To obtain sub-second accuracy in sidereal time, it is important
# to use the form of Universal Time called UT1 as the basis for the input
# Julian date.
# The maximum value of the equation of the equinoxes is about 1.1 seconds, so
# if an error of ~1 second is unimportant, the last series of formulas can
# be skipped entirely.  In this case set eqeq = 0 and GAST = GMST, and use
# either UT1 or UTC as the Universal Time basis for the input Julian date.
#
