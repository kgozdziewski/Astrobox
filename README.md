# Astrobox

A very simple (almost trivial) Python 3.x package for introductory Spherical and Classical Astronomy course. It focuses on checking trigonometric formulae developed in scalar form. Note that some functions may lack (sometimes important) details, be incomplete or oversimplified. This version has a Polish interface only.

Bardzo prosty (trywialny) pakiet Pythona 3.x dla podstawowego kursu Astronomii Sferycznej i Klasycznej. Skupia się na wzorach trygonometrycznych opracowanych w postaci skalarnej. Niektórym funkcjom może brakować (czasem ważnych) szczegółów, mogą one być niekompletne lub nadmiernie uproszczone do poważniejszych zastosowań. 

Astrobox korzysta ze standardowego modułu do obliczeń matematycznych numpy
  
    * Użycie: 
      import astrobox as asa
      (funkcje modułu należy następnie poprzedzić prefiksem asa)
  
    * Lista funkcji w module: 
      dir(asa)
  
    * Opis funkcji "foo" z przykładem: 
      asa.opis(asa.foo)
      asa.doc(asa.foo)
  
    * Uaktualnienie modułu po zmianach w sesji Python'a:
      from importlib import reload
      reload(asa)
      
© Krzysztof Goździewski, 2009-2023, wersja 18.03.2023
