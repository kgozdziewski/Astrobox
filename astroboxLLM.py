# -*- coding: utf-8 -*-
"""
astrobox - moduł do obliczeń astronomii sferycznej i klasycznej

Użytkowanie:
    import astrobox as asa
    asa.<funkcja>(...)      # funkcje poprzedzamy prefiksem asa
    asa.test()              # uruchomienie testów
    asa.doc(asa.<funkcja>)  # opis funkcji

Aktualizacja w powłoce Pythona po zmianach w pliku:
    from importlib import reload
    reload(asa)

© Krzysztof Goździewski, 2009-2026, wersja 24.03.2026
"""

import numpy as np
import inspect

# =============================================================================
# Stałe fundamentalne
# =============================================================================
eps    = 2.26e-16       # zero maszynowe (dokładność 1+eps=1)
d2r    = np.pi / 180.0  # czynnik konwersji stopni na radiany
JD2000 = 2451545.0      # data juliańska epoki J2000.0 (January 1.5, 2000 TT)
cy     = 36525.0        # stulecie juliańskie

# =============================================================================
# Konwersje kątów (radiany <-> czasowe/stopniowe)
# =============================================================================

def hr2rad(h, m, s):
    """
    Zamiana kąta w mierze czasowej na łukową (rad).
    
    Parametry:
        h, m, s — godziny, minuty, sekundy
    
    Zwraca:
        kąt w radianach
    
    Przykład:
        >>> round(hr2rad(20, 30, 30), 9)
        5.369069111
    """
    return (h + m/60.0 + s/3600.0) * np.pi / 12.0

def deg2rad(d, m, s):
    """
    Zamiana kąta w mierze stopniowej na łukową (rad).
    
    Parametry:
        d, m, s — stopnie, minuty, sekundy
    
    Zwraca:
        kąt w radianach
    
    Przykład:
        >>> round(deg2rad(20, 30, 30), 9)
        0.357937941
    """
    return (d + m/60.0 + s/3600.0) * np.pi / 180.0

def _angle_to_dms(rad, to_deg=True):
    """
    Wewnętrzna funkcja: zamiana radianów na (stopnie/godziny, minuty, sekundy).
    
    Parametry:
        rad — kąt w radianach
        to_deg — True dla stopni (180°=π), False dla godzin (12h=π)
    """
    factor = 180.0/np.pi if to_deg else 12.0/np.pi
    angle = np.mod(rad, 2*np.pi) * factor
    d, deg = np.modf(angle)
    m, sec = np.modf(d * 60.0)
    return deg, m, sec * 60.0

def rad2deg(rad):
    """
    Zamiana kąta w mierze łukowej na stopniową (deg, min, sec).
    
    Zwraca:
        tuple (stopnie, minuty, sekundy)
    
    Przykład:
        >>> rad2deg(0.123456)
        (7.0, 4.0, 24.627920047527176)
    """
    return _angle_to_dms(rad, to_deg=True)

def rad2hr(rad):
    """
    Zamiana kąta w mierze łukowej na czasową (godziny, minuty, sekundy).
    
    Zwraca:
        tuple (godziny, minuty, sekundy)
    
    Przykład:
        >>> rad2hr(0.123456)
        (0.0, 28.0, 17.64186133610181)
    """
    return _angle_to_dms(rad, to_deg=False)

def day2hr(dt):
    """
    Zamiana interwału czasu na (dzień, godziny, minuty).
    
    Parametry:
        dt — liczba dni (może być ułamkowa)
    
    Zwraca:
        tuple (dzień, godzina, minuty)
    
    Przykład:
        >>> day2hr(1.5)
        (1, 12.0, 0.0)
    """
    day = int(dt)
    h, hr = np.modf((dt - day) * 24)
    return day, hr, h * 60.0

# =============================================================================
# Formatowanie kątów — wypisywanie w formacie DMS
# =============================================================================

def _format_angle(rad, unit, name='', to_deg=True, fsec='%3.2f'):
    """
      funkcja formatująca kąt (ułatwiająca powtórzenia).
    """
    
    d, m, sec = _angle_to_dms(rad, to_deg=to_deg)
    
    # normalizacja do [-180°, 180°] lub [0°, 360°] zgodnie z kontekstem
    if to_deg and abs(rad) > np.pi:
        d, m, sec = _angle_to_dms(rad - 2*np.pi, to_deg=True)
    
    # format wyjściowy
    if unit == 'deg':
        suffix = "° %2.0f'" + fsec + "\""
    elif unit == 'hr':
        suffix = "h %2.0fm " + fsec + "s"
    else:
        raise ValueError("unit must be 'deg' or 'hr'")
    
    if name:
        print(f"{name} = {d:.0f}{suffix}" % (m, sec))
    else:
        print(f"{d:.0f}{suffix}" % (m, sec))

def printdeg(angle, name='', fsec='%3.2f'):
    """
    wypisuje kąt w stopniach (DMS).
    
    parametry:
        angle - kąt w radianach
        name - opcjonalna nazwa zmiennej
        fsec - format sekund
    
    przykład:
        >>> printdeg(np.pi/2, 'pół obrotu')
        pół obrotu = 90°  0' 0.00"
    """
    # przyjmujemy kąt modulo 360°
    angle = np.mod(angle, 2*np.pi)
    d, m, sec = _angle_to_dms(angle, to_deg=True)
    
    if name:
        print(f"{name} = {d:.0f}° {m:.0f}'{fsec}\"")
    else:
        print(f"{d:.0f}° {m:.0f}'{fsec}\"")

def printhr(angle, name='', fsec='%3.2f'):
    """
    wypisuje kąt w godzinach (HMS).
    
    parametry:
        angle - kąt w radianach
        name - opcjonalna nazwa zmiennej
        fsec - format sekund
    
    Przykład:
        >>> printhr(np.pi, 'pół obrotu')
        pół obrotu = 12h  0m  0.00s
    """
    angle = np.mod(angle, 2*np.pi)
    h, m, sec = _angle_to_dms(angle, to_deg=False)
    
    if name:
        print(f"{name} = {h:.0f}h {m:.0f}m {fsec}s")
    else:
        print(f"{h:.0f}h {m:.0f}m {fsec}s")

def printfhr(dt, name=''):
    """
    wypisuje datę juliańską w formacie kalendarzowym (do minut).
    
    parametry:
        dt - data juliańska
        name - opcjonalna nazwa
    
    przykład:
        >>> printfhr(2450000.1, 'test')
        test = 1995 10  9 14h 24.0m
    """
    y, mo, day = JDtoDate(dt)
    d, hr, m = day2hr(day)
    if name:
        print(f"{name} = {y:.0f} {mo:.0f} {d:.0f} {hr:.0f}h {m:.1f}m")
    else:
        print(f"{y:.0f} {mo:.0f} {d:.0f} {hr:.0f}h {m:.1f}m")

# =============================================================================
# Współrzędne: równikowe ↔ horyzontalne
# =============================================================================

def eqI2hor(phi, t, delta):
    """
    przemiana współrzędnych równikowych I (α, δ) na horyzontalne (A, h).
    
    parametry:
        phi - szerokość geograficzna obserwatora (rad)
        t - kąt godzinny (rad)
        delta - deklinacja (rad)
    
    zwraca:
        (A, h) — azymut, wysokość (rad)
    
    uwaga:
        funkcja nie działa poprawnie dla obserwatora zbyt blisko bieguna.
    """
    if abs(abs(phi) - np.pi/2) < eps:
        print("# <eqI2hor>: obserwator zbyt blisko bieguna!")
        return 0.0, 0.0
    
    sh = np.sin(delta)*np.sin(phi) + np.cos(delta)*np.cos(phi)*np.cos(t)
    if abs(abs(sh) - 1) < eps:
        print("# <eqI2hor>: wysokość zbyt bliska ±90°!")
        return 0.0, 0.0
    
    h  = np.arcsin(sh)
    cA = (np.sin(delta) - np.sin(h)*np.sin(phi)) / (np.cos(h)*np.cos(phi))
    sA = np.sin(t)*np.cos(delta) / np.cos(h)
    A  = np.arctan2(sA, cA)
    return A, h

def hor2eqI(phi, A, h):
    """
    przemiana współrzędnych horyzontalnych (A, h) na równikowe I (t, δ).
    
    parametry:
        phi - szerokość geograficzna (rad)
        A - azymut (rad)
        h - wysokość (rad)
    
    zwraca:
        (t, delta) — kąt godzinny, deklinacja (rad)
    """
    if abs(abs(phi) - np.pi/2) < eps:
        print("# <hor2eqI>: obserwator zbyt blisko bieguna!")
        return 0.0, 0.0
    
    sd = np.sin(h)*np.sin(phi) + np.cos(h)*np.cos(phi)*np.cos(A)
    if abs(abs(sd) - 1) < eps:
        print("# <hor2eqI>: deklinacja zbyt bliska ±90°!")
        return 0.0, 0.0
    
    delta = np.arcsin(sd)
    ct = (np.sin(h) - np.sin(phi)*np.sin(delta)) / (np.cos(phi)*np.cos(delta))
    st = np.sin(A)*np.cos(h) / np.cos(delta)
    t = np.arctan2(st, ct)
    return t, delta

def eqII2hor(phi, LST, alpha, delta):
    """
    przemiana współrzędnych równikowych II (α, δ) na horyzontalne (A, h).
    
    parametry:
        phi - szerokość geograficzna (rad)
        LST - lokalny czas gwiazdowy (rad)
        alpha, delta - współrzędne równikowe II (rad)
    
    zwraca:
        (A, h) - azymut, wysokość (rad)
    
    uwaga:
        używa przekształcenia t = LST - α, potem eqI2hor.
    """
    t = LST - alpha
    return eqI2hor(phi, t, delta)

def hor2eqII(phi, LST, A, h):
    """
    przemiana współrzędnych horyzontalnych (A, h) na równikowe II (α, δ).
    
    parametry:
        phi, LST, A, h — jak wyżej
    
    zwraca:
        (alpha, delta) — współrzędne równikowe II (rad)
    """
    t, delta = hor2eqI(phi, A, h)
    return LST - t, delta

# =============================================================================
# Współrzędne: równikowe II <-> ekliptyczne
# =============================================================================

def ecl2eqII(epsilon, lam, beta):
    """
    przemiana współrzędnych ekliptycznych (λ, β) na równikowe II (α, δ).
    
    parametry:
        epsilon - nachylenie ekliptyki do równika (rad)
        lam, beta - długość i szerokość ekliptyczna (rad)
    
    zwraca:
        (alpha, delta) - rektascensja i deklinacja (rad)
    """
    sd = np.sin(beta)*np.cos(epsilon) + np.cos(beta)*np.sin(epsilon)*np.sin(lam)
    if abs(abs(sd) - 1) < eps:
        print("# <ecl2eqII>: deklinacja zbyt bliska ±90°!")
        return 0.0, 0.0
    
    delta = np.arcsin(sd)
    sa = (-np.sin(beta) + np.sin(delta)*np.cos(epsilon))/(np.cos(delta)*np.sin(epsilon))
    ca = np.cos(lam)*np.cos(beta) / np.cos(delta)
    alpha = np.arctan2(sa, ca)
    return alpha, delta

def eqII2ecl(epsilon, alpha, delta):
    """
    przemiana współrzędnych równikowych II (α, δ) na ekliptyczne (λ, β).
    
    parametry:
        epsilon - nachylenie ekliptyki do równika (rad)
        alpha, delta - rektascensja i deklinacja (rad)
    
    zwraca:
        (lam, beta) - długość i szerokość ekliptyczna (rad)
    """
    sb = np.sin(delta)*np.cos(epsilon) - \
         np.cos(delta)*np.sin(epsilon)*np.sin(alpha)
    if abs(abs(sb) - 1) < eps:
        print("# <eqII2ecl>: szerokość ekliptyczna zbyt bliska ±90°!")
        return 0.0, 0.0
    
    beta = np.arcsin(sb)
    sl = (np.sin(delta) - np.sin(beta)*np.cos(epsilon)) / (np.cos(beta)*np.sin(epsilon))
    cl = np.cos(alpha)*np.cos(delta) / np.cos(beta)
    lam = np.arctan2(sl, cl)
    return lam, beta

# =============================================================================
# Daty juliańskie i czas gwiazdowy
# =============================================================================

def JDN(year, month, day):
    """
    dzień juliański (JDN) - całkowita część JD.
    
    parametry:
        year, month, day - data kalendarzowa (UTC)
    
    zwraca:
        JDN jako float
    
    uwaga:
        działa poprawnie dla dat od 15 października 1582 (kalendarz gregoriański).
    
    przykład:
        >>> JDN(2000, 1, 1)
        2451545.0
    """
    a = np.floor((14 - month)/12.)
    y = year + 4800 - a
    m = month + 12*a - 3
    return day + np.floor((153*m + 2)/5.) + y*365 + \
           np.floor(y/4.) - np.floor(y/100.) + np.floor(y/400.) - 32045.0

def JD(year, month, day, hour, minute, second):
    """
    data juliańska (JD).
    
    parametry:
        year, month, day, hour, minute, second — kalendarzowa data i czas
    
    zwraca:
        JD jako float
    
    przykład:
        >>> JD(2000, 1, 1, 12, 0, 0)
        2451545.0
    """
    return JDN(year, month, day) + (hour - 12.0)/24.0 + minute/1440.0 + second/86400.0

def GMST(year, month, day, hour, minute, second):
    """
    średni czas gwiazdowy Greenwich (GMST) w radianach.
    
    parametry:
        year, month, day, hour, minute, second — data UT1/UTC
    
    zwraca:
        GMST (rad), dokładność ~1s w latach 1900–2100
    
    przykład:
        >>> round(rad2hr(GMST(2010,6,10,17,30,30))[2], 3)
        20.459
    """
    D = JD(year, month, day, hour, minute, second) - JD2000
    time = np.mod(18.697374558 + 24.06570982441908*D, 24)
    return time * 15 * np.pi / 180.0

def GST(year, month, day, hour, minute, second):
    """
    prawdziwy czas gwiazdowy Greenwich (GAST) w radianach.
    
    parametry:
        jak wyżej
    
    zwraca:
        GAST (rad), uwzględnia nutację i równanie równicy
    
    przykład:
        >>> round(rad2hr(GST(2019,3,13,0,0,0))[2], 2)
        21.02
    """
    D = JD(year, month, day, hour, minute, second) - JD2000
    time = np.mod(18.697374558 + 24.06570982441908*D, 24)
    
    # Nutacja i równanie równicy
    epsilon = np.pi * (23.4393 - 0.0000004*D) / 180.0
    Omega = np.pi * (125.04 - 0.052954*D) / 180.0
    L = np.pi * (280.47 + 0.98565*D) / 180.0
    dPsi = -0.000319*np.sin(Omega) - 0.000024*np.sin(2*L)
    eqeq = dPsi * np.cos(epsilon)
    
    return (time + eqeq) * 15 * np.pi / 180.0

def JDtoDate(jd):
    """
    xamiana daty Juliańskiej na kalendarzową.
    
    parametry:
        jd - data Juliańska
    
    xwraca:
        (year, month, day) - rok, miesiąc, dzień (float, może być ułamkowy)
    
    przykład:
        >>> JDtoDate(2451545.0)
        (2000.0, 1.0, 1.5)
    """
    jd += 0.5
    F, I = np.modf(jd)
    I = int(I)
    A = np.trunc((I - 1867216.25)/36524.25)
    B = I + 1 + A - np.trunc(A/4.) if I > 2299160 else I
    C = B + 1524
    D = np.trunc((C - 122.1)/365.25)
    E = np.trunc(365.25 * D)
    G = np.trunc((C - E)/30.6001)
    day = C - E + F - np.trunc(30.6001 * G)
    
    month = G - 1 if G < 13.5 else G - 13
    year = D - 4716 if month > 2.5 else D - 4715
    return year, month, day

# =============================================================================
# konwersje współrzędnych wektorowych: biegunowe <-> kartezjańskie
# =============================================================================

def sph2cart(r, lon, phi):
    """
    lonwersja biegunowych [r, lon, phi] do kartezjańskich [x, y, z].
    
    parametry:
        r - promień
        lon - długość (rad)
        phi - szerokość (rad)
    
    zwraca:
        (x, y, z) - współrzędne kartezjańskie
    
    przykład:
        >>> sph2cart(1, 0, 0)
        (1.0, 0.0, 0.0)
    """
    return r * np.cos(phi) * np.cos(lon), \
           r * np.cos(phi) * np.sin(lon), \
           r * np.sin(phi)

def cart2sph(x, y, z):
    """
    konwersja współrzędnych kartezjańskich [x, y, z] 
    do biegunowych [r, lon, phi].
    
    parametry:
        x, y, z - współrzędne kartezjańskie
    
    zwraca:
        (r, lon, phi)
    """
    rho = np.sqrt(x*x + y*y + z*z)
    lon = np.mod(np.arctan2(y/rho, x/rho), 2*np.pi)
    phi = np.mod(np.arcsin(z/rho), 2*np.pi)
    return rho, lon, phi

def sph2cartv(r, lon, phi):
    """
      Wersja wektorowa sph2cart - zwraca numpy array.
    """
    return np.array([r * np.cos(phi) * np.cos(lon),
                     r * np.cos(phi) * np.sin(lon),
                     r * np.sin(phi)])

def cart2sphv(v):
    """
      wersja wektorowa cart2sph — przyjmuje numpy array.
    """
    rho = np.linalg.norm(v)
    lon = np.mod(np.arctan2(v[1]/rho, v[0]/rho), 2*np.pi)
    phi = np.mod(np.arcsin(v[2]/rho), 2*np.pi)
    return rho, lon, phi

# =============================================================================
# macierze obrotu
# =============================================================================

def A1(phi):
    """
      macierz obrotu wokół osi x (1) o kąt phi.
    """
    return np.array([[1, 0, 0],
                     [0, np.cos(phi), np.sin(phi)],
                     [0, -np.sin(phi), np.cos(phi)]])

def A2(theta):
    """
      macierz obrotu wokół osi y (2) o kąt theta.
    """
    return np.array([[np.cos(theta), 0, -np.sin(theta)],
                     [0, 1, 0],
                     [np.sin(theta), 0, np.cos(theta)]])

def A3(psi):
    """
      macierz obrotu wokół osi z (3) o kąt psi.
    """
    return np.array([[np.cos(psi), np.sin(psi), 0],
                     [-np.sin(psi), np.cos(psi), 0],
                     [0, 0, 1]])

def precess(dt):
    """
    macierz precesji ogólnej (IAU 1986).
    
    parametry:
        dt - różnica czasu w dniach od J2000.0
    
    zwraca:
        macierz 3x3 - przekształca wektor z epoki t0 na t0+dt
    """
    t = dt / cy
    zA = (0.6406161*t + 0.0003041*t**2 + 5e-7*t**3) * d2r
    thetaA = (0.5567530*t - 1.18e-5*t**2 - 1.1e-6*t**3) * d2r
    ksiA = (0.6406161*t + 8.39e-5*t**2 + 5e-6*t**3) * d2r
    return A3(-zA) @ A2(thetaA) @ A3(-ksiA)

# =============================================================================
# współrzędne Słońca i równanie czasu
# =============================================================================

def SunLongitude(year, month, day, hour, min_, sec):
    """
    długość ekliptyczna Słońca (λ☉) w radianach.
    
    algorytm: Almagest, Fitzpatrick (dokładność ~1').
    uwzględnia: precesję punktu Barana, precesję peryhelium, aberrację.
    
    parametry:
        jak JD
    
    zwraca:
        λ☉ (rad)
    """
    t = JD(year, month, day, hour, min_, sec)
    ecc = 0.01671123
    omega0 = 102.93768193
    kappa = 20.5/3600.0
    psi = 3.8246e-5
    nomega = 0.32327364/36525.0
    n = 35999.37244981/36525.0
    lambda0 = 100.46457166 - kappa
    
    meanL = np.mod(lambda0 + (n + psi)*(t - JD2000), 360.0)
    M = np.mod(lambda0 - omega0 + (n - nomega)*(t - JD2000), 360.0)
    q = 2*ecc*np.sin(M*d2r) + 1.25*ecc**2*np.sin(2*M*d2r)
    
    return np.mod(meanL + q + 180.0, 360.0) * d2r

def SunLambda(year, month, day, hour, min_, sec):
    """
      długość ekliptyczna Słońca (uproszczona, Astronomical Almanac)
      dokładność: ~1' w latach 1800–2050
    """
    dt = JD(year, month, day, hour, min_, sec) - JD2000
    meanL = np.mod(280.460 + 0.9856474*dt, 360.0)
    M = np.mod(357.528 + 0.9856003*dt, 360.0)
    q = 1.915*np.sin(M*d2r) + 0.020*np.sin(2*M*d2r)
    return np.mod(meanL + q, 360.0) * d2r

def SunLongitudet(jd):
    """
      SunLongitude z argumentem JD (użyteczne dla tablic).
    """
    ecc = 0.01671123
    omega0 = 102.93768193
    kappa = 20.5/3600.0
    psi = 3.8246e-5
    nomega = 0.32327364/36525.0
    n = 35999.37244981/36525.0
    lambda0 = 100.46457166 - kappa
    
    meanL = np.mod(lambda0 + (n + psi)*(jd - JD2000), 360.0)
    M = np.mod(lambda0 - omega0 + (n - nomega)*(jd - JD2000), 360.0)
    q = 2*ecc*np.sin(M*d2r) + 1.25*ecc**2*np.sin(2*M*d2r)
    
    return np.mod(meanL + q + 180.0, 360.0) * d2r

def DeltaE(year, month, day, hour, min_, sec):
    """
      równanie czasu - różnica między czasem prawdziwym a średnim (minuty).
    
    parametry:
        jak JD
    
    zwraca:
        ΔE w minutach (T_true - T_mean)
    """
    lsun = SunLongitude(year, month, day, hour, min_, sec)
    M = np.mod(lsun + (357.588 - 280.458)*d2r, 2*np.pi)
    eps = deg2rad(23, 26, 21.45)
    
    # rektascensja Słońca
    asun = np.arctan(np.cos(eps)*np.tan(lsun))
    if lsun >= 3*np.pi/2:
        asun += 2*np.pi
    elif np.pi/2 < lsun < 3*np.pi/2:
        asun += np.pi
    
    dalpha = lsun - asun - 2*0.01671123*np.sin(M)
    return dalpha * 180/np.pi / 15 * 60

# =============================================================================
# refrakcja atmosferyczna
# =============================================================================

def refract(zobs, Hobs=0, P=1013, T=0):
    """
      oblicza refrakcję atmosferyczną w kącie zenitalnym (rad).
    
     algorytm: NOVAS (USNO/AA), Bennett (1982).
    
    parametry:
        zobs - odległość zenitalna obserwowana (rad), musi być w [0.1, 91°]
        Hobs - wysokość obserwatora nad morzem (m)
        P - ciśnienie atmosferyczne (hPa)
        T - temperatura (°C)
    
    zwraca:
        refrakcja (rad)
    
    uwaga:
        nie nadaje się do precyzyjnej redukcji obserwacji.
    """
    if zobs < 0.1 or zobs > 91*np.pi/180:
        return 0.0
    
    h = 90.0 - zobs*180/np.pi  # wysokość w stopniach
    
    # ciśnienie rzeczywiste (skalowanie z wysokości)
    p = P * np.exp(Hobs/9100) if P == 1010 and T == 0 else P
    
    # wzór Bennett
    r = 0.016667 / np.tan((h + 7.31/(h + 4.4)) * np.pi/180.)
    return r * (0.28 * p / (T + 273.0)) * np.pi/180.
    
def MoonLongitude(jd):
    """
    długość i szerokość ekliptyczna Księżyca (ELP-2000, wersja uproszczona).
    
    dokładność:
        λ: ~10–15", β: ~5–7" (do przybliżonych obliczeń).
    
    parametry:
        jd - data Juliańska (float)
    
    zwraca:
        (lambda_m, beta_m) - długość i szerokość ekliptyczna (rad)
    
    uwaga:
        model zawiera główne czynniki: ruch średni, librację, precesję, nutację
        oraz główne okresy (anomalia, węzeł, perigeum).
    """
    # interwał od J2000.0 w dniach
    T = (jd - JD2000) / 36525.0  # wieki juliańskie
    
    # ruch średni (stopnie/dzień)
    n = 13.1763966  
    # średnie ruchy: M_L (długość), M_S (anomalia), F (argument szerokości)
    # stałe (stopnie)
    L0 = 270.434164    # średnia długość
    M0 = 358.616448    # średnia anomalies
    C0 = 299.770617    # średnia długość węzła
    P0 = 345.466868    # średnia długość perigeum
    
    # kąty w stopniach
    L = np.mod(L0 + 13.1763966 * (jd - JD2000), 360.0)
    M = np.mod(M0 + 0.1113957 * (jd - JD2000), 360.0)
    F = np.mod(C0 + 13.0649930 * (jd - JD2000), 360.0)
    omega = np.mod(P0 + 0.0034730 * (jd - JD2000), 360.0)
    
    # zamiana na radiany
    L, M, F, omega = np.radians([L, M, F, omega])
    
    # Głównie wyrazowy ruch w długości (przybliżenie ELP-2000)
    # wyraz 1: 6.289°(sin M)
    # wyraz 2: 0.214°(sin 2M)
    # wyraz 3: 0.658°(sin 2(L - F)) — wariacja
    # wyraz 4: 0.185°(sin M') — M' = 2(L - F) - M
    # wyraz 5: 0.114°(sin 2F)
    # wyraz 6: 0.041°(sin M'')
    M_prime = 2*(L - F) - M
    M_second = 2*(L - F)
    
    delta_L = (6.289*np.sin(M) + 
               0.214*np.sin(2*M) + 
               0.658*np.sin(2*(L - F)) + 
               0.185*np.sin(M_prime) + 
               0.114*np.sin(2*F) + 
               0.041*np.sin(M_prime + M))
    
    # Średnia długość + poprawki
    lambda_m = np.radians(np.mod(L0 + 13.1763966*(jd-JD2000) + delta_L, 360.0))
    
    # Szerokość ekliptyczna (w przybliżeniu)
    #主要 wyrazy:
    # 5.128°sin F
    # 0.281°sin(F + M)
    # 0.278°sin(F - M)
    # 0.173°sin(F - 2(L - F))
    # 0.055°sin(F - 2M)
    beta_m_rad = (5.128*np.sin(F) + 
                  0.281*np.sin(F + M) + 
                  0.278*np.sin(F - M) + 
                  0.173*np.sin(F - 2*(L - F)) + 
                  0.055*np.sin(F - 2*M))
    beta_m = np.radians(beta_m_rad)
    
    return lambda_m, beta_m    
    
def MoonAltitude(phi, lambda_loc_deg, jd):
    """
    oblicza wysokość (i azymut) Księżyca w układzie horyzontalnym.
    
    parametry:
        phi - szerokość geograficzna obserwatora (rad)
        lambda_loc_deg - długość geograficzna obserwatora (°, dodatnia = wschód)
        jd - data Juliańska
    
    zwraca:
        (A, h) — azymut, wysokość (rad)
    
    uwaga:
        współrzędne Księżyca obliczane są w epoce JD, 
        z uwzględnieniem zmiany nachylenia ekliptyki.
    """
    # 1. nachylenie ekliptyki w epoce JD (IAU 1980)
    T = (jd - JD2000) / cy
    epsilon = np.radians(23.43929111 - 0.013004167*T - 1.64e-7*T**2 + 5.04e-7*T**3)
    
    # 2. pozycja Księżyca w układzie ekliptycznym
    lambda_m, beta_m = MoonLongitude(jd)
    
    # 3. przelicz (λ, β) → (α, δ)
    alpha_m, delta_m = ecl2eqII(epsilon, lambda_m, beta_m)
    
    # 4. oblicz czas gwiazdowy (LST)
    # najpierw JD → rok, miesiąc, dzień, godzina
    y, mo, day = JDtoDate(jd)
    d, hr, m = day2hr(day)
    sec = m * 60
    gmst = GMST(int(y), int(mo), int(day), int(hr), int(m), sec)
    
    # LST = GMST + długość lokalna (w godzinach: λ/15)
    # λ w stopniach → radiany (1h=15°)
    lst = gmst + np.radians(lambda_loc_deg / 15.0 * 15.0)  
    
    # 5. przelicz na współrzędne horyzontalne
    A, h = eqII2hor(phi, lst, alpha_m, delta_m)
    
    return A, h    

# =============================================================================
# Funkcje pomocnicze dla studentów
# =============================================================================

def opis(fun):
    """Krótki opis funkcji."""
    print(f"\nOpis funkcji {fun}")
    print(inspect.getdoc(fun))

def doc(fun):
    """
      dokumentacja funkcji (to samo co opis).
    """
    print(inspect.getdoc(fun))

# =============================================================================
# Testy
# =============================================================================

def test():
    """
      Przykładowe testy modułu.
    """
    print("\n=== Przykładowe testy astrobox ===\n")
    
    # deklinacja Polaris
    DEC = deg2rad(89, 15, 40.8)
    printdeg(DEC, 'deklinacja Polaris')
    
    # czas gwiazdowy
    gmst = GMST(2019, 3, 18, 14, 0, 0)
    printhr(gmst, 'GMST')
    gst = GST(2019, 3, 18, 14, 0, 0)
    printhr(gst, 'GAST')
    
    # przemiana ekliptyczna -> równikowa
    alpha, delta = ecl2eqII(23.5*d2r, 30.0*d2r, 23.5*d2r)
    print(f"\nα = {alpha:.3f} rad, δ = {delta:.3f} rad")
    
    # Słońce w równonocy
    print("\nDługość Słońca w równonocy 20.III (Almagest):")
    for y, m, d, h, mi in [
        (1896, 3, 20, 2, 46), (2010, 3, 20, 17, 32),
        (2015, 3, 20, 22, 45), (2020, 3, 20, 3, 49)]:
        printdeg(SunLongitude(y, m, d, h, mi, 0), f"λ☉ {y}")
    
    # równanie czasu (tygodnie 2026)
    print("\nRównanie czasu w 2026 (w minutach):")
    year = 2026
    days_in_month = [31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31]
    for mo in range(1, 13):
        for d in range(1, days_in_month[mo-1], 7):
            de = DeltaE(year, mo, d, 12, 0, 0)
            print(f"{mo:2d}/{d:2d}: {de:6.2f} min")
            
    # Księżyc
    print( "\nwspółrzędme Księzyca w układzie ekliptycznym")
    jd = JD(2024, 6, 4, 12, 0, 0)
    lam, beta = MoonLongitude(jd)
    printdeg(lam, "λ Moon")
    printdeg(beta, "β Moon")
    pos = MoonPosition(jd)
    print(f"wpółrzędne Księżyca (AU): {pos}")            
    
    # Obserwator w Warszawie (52°07'N, 20°59'E)
    lat = 52 + 7/60.0
    lon = 20 + 59/60.0

    # Data: 4 czerwca 2024, 22:00 UT
    A, h = MoonAltitudeSimple(lat, lon, 2024, 6, 4, 22, 0, 0)

    print("wysokość Księżyca:")
    printdeg(h, "h")
    print("azymut:")
    printdeg(A, "A")    

# ==============================================================================
# Uruchomienie przy starcie modułu (opcjonalne)
# ==============================================================================
if __name__ == "__main__":
    print(f"Moduł astrobox © Krzysztof Goździewski {__copyright__}, wersja {__version__}")
    import doctest
    print("\ntesty:")
    doctest.testmod(optionflags=doctest.ELLIPSIS)
    test()
