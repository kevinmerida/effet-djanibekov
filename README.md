# L'effet Djanibekov

## Description

L'effet se produit pour tout corps rigide en impesanteur — ou chute libre — qui présente trois axes principaux d'inerties différentes, selon la [description](https://fr.wikipedia.org/wiki/Effet_Djanibekov) disponible sur Wikipédia.

Ce phénomène est illustré dans une [vidéo](https://www.youtube.com/watch?v=SAQ-iIJkLzA&t=277s) sur la chaîne **Véritasium** (le lien inclue le time code qui mène directement à la séquence montrant le phénomène).

Une [première vidéo](https://www.youtube.com/watch?v=BzJsEE4yTJw) sur la chaîne **Kevin Merida** présentait une simulation numérique du phénomène.


Pour améliorer notablement la précision sur la simulation numérique, la méthode décrite dans le [document](https://www.f-legrand.fr/scidoc/srcdoc/sciphys/meca/solide/solide-pdf.pdf) rédigé par [Frédéric Legrand](https://www.f-legrand.fr/scidoc/) a été utilisée.

Les résultats de cette nouvelle simulation et les programmes associés sont accessibles dans ce notebook.

## Mise en équation

On définit le vecteur rotation $\vec{\Omega}_s$ dans le repère solidaire du corps rigide. Ce repère correspond aux axes principaux d'inertie, et les moments principaux d'inertie y sont notés $I_1$, $I_2$ et $I_3$.

$$
\vec{\Omega}_s=\left(
\begin{matrix}
\Omega_1\\
\Omega_2\\
\Omega_3
\end{matrix}\right)
$$

Par ailleurs, on définit les trois vecteurs orthonormaux orientant les axes principaux d'inertie du corps rigide, rassemblés dans les trois colonnes d'une matrice $Q$. Le mouvement s'effectue dans un référentiel inertiel et aucune force ne s'applique au solide.

L'équation dynamique est donnée ci-dessous. Elle comprend 12 variables, 3 pour les vitesses angulaires et 9 pour la matrice Q qui permettra de représenter le mouvement du solide dans le référentiel inertiel.

$$
\begin{align}
\frac{d\Omega_1}{dt}&=\Omega_2\Omega_3\frac{(I_2 - I_3)}{I_1}\\
\frac{d\Omega_2}{dt}&=\Omega_3\Omega_1\frac{(I_3 - I_1)}{I_2}\\
\frac{d\Omega_3}{dt}&=\Omega_1\Omega_2\frac{(I_1 - I_2)}{I_3}\\
\frac{dQ}{dt}&=Q\left(\begin{matrix}
0&-\Omega_3&\Omega_2\\
\Omega_3&0&-\Omega_1\\
-\Omega_2&\Omega_1&0
\end{matrix}\right)\\
\end{align}
$$

La résolution numérique de cette équation dynamique peut être validée en calculant le moment cinétique $\vec{\sigma}$ dans le référentiel inertiel, qui doit rester constant en théorie. Le moment cinétique $\vec{\sigma}_s$ défini sur les axes principaux d'inertie a une expression simple et donne accès à $\vec{\sigma}$ via le changement de base réalisé avec la matrice $Q$.

$$
\vec{\sigma}=Q\vec{\sigma}_s=Q\left(
\begin{matrix}
I_1 \Omega_1\\
I_2\Omega_2\\
I_3\Omega_3
\end{matrix}
\right)
$$

## Programmes et résultats

Le programme [**toupie.py**](Code/toupie.py) contient les fonctions permettant de réaliser les simulations numériques, puis de visualiser les résultats. Le notebook [**toupie.ipynb**](Notebook/toupie.ipynb) donne quelques exemples d'appels à ces fonctions.

https://github.com/user-attachments/assets/4b143d46-965e-40d2-8ee5-ca98e5521ed6
