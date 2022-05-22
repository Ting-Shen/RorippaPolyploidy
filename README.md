This is the R scripts used in the project

asr_ploidy.R

The R script used for ancestral ploidy state reconstruction under two or multiple ploidy level models (TPL or MPL model).

In the custom TPL model, species were divided into diploidy or polyploidy according to the basic chromosome number x = 8. Transition can only happen from diploidy (state x) to polyploidy (state y) but not the reverse:

![image](https://user-images.githubusercontent.com/48637894/169685702-59f9a106-f34b-4758-9971-fc8c3210c5b2.png)

In the custom MPL model, species were divided into multiple ploidy levels ranging from diploidy (2x; x means the basic chromosome number) to high ploidy (-x; the number (-) before x means the corresponding ploidy levels). The chromosome number information for Rorippa species was obtained using information from available databases or publications described above. Transition can only happen from low (state group of x) to high ploidy (state group of y) 
under equal rate model (MPL-ER):

![image](https://user-images.githubusercontent.com/48637894/169685589-bec51aac-3a29-428d-bb15-4742ba223410.png)

or under all rates different model (MPL-ARD):

![image](https://user-images.githubusercontent.com/48637894/169685627-eab838f2-2e82-40b6-ab39-ac8c2f94dc6f.png)

or under symmetric model (MPL-SYM):

![image](https://user-images.githubusercontent.com/48637894/169685638-67863d13-d0e2-48c6-83bd-05e0e0fb4023.png)
