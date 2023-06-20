# BoneEstrogenHPA

A combined mathematical model of the hypothalamus-pituitary-adrenal (HPA) axis and the calcium-
parathyroid-vitamin D (Ca-PTH-D) axis with added effects of mean estradiol in premenopausal women, decreased estradiol according to postmenopausal levels, and ethinyl estradiol according to levels in combined oral contraceptives.

The model was made as part of a master thesis project where it was simulated with a stress input, under primary hyperparathyroidism, under secondary hyperparathyroidism or under primary hypoparathyroidism.

The Bone model by Peterson and Riggs from 2010 [[1]](https://doi.org/10.1016/j.bone.2009.08.053) serves as the foundation, it is expanded with the incorporation of the HPA axis model developed by Karin et al. from 2020 [[2]](https://doi.org/10.15252/msb.20209510).

## How to run it
Simulations are created by running the main.r file in the model_code folder. Plots of the simulations can then be rendered by running the create_simulation_plots rmd. 

## How to change it
Constant parameters can be added to params.r. 
Introduction of new species to the model must have an initial value in the initial.r script as well as a name in the compartments.r script.
Differential equations for the new species are added in the model.r script together with other effects.
After making changes to model.r, a new steady state of this new model must me created before simulating stress or parathyroid disease states. This will probably take a while, depending on computer performance.


## Introduced effects
- Stimulatory effect of ethinyl estradiol (EE) on _total serum_ cortisol
- Stimulatory effect of PTH on cortisol (no clear evidence for this)
- Effects of estrogen (independently of type) on different bone factors

## Limitations
There are a lot of limitations to this model, some of them are:
- It does not contain any age-related effects.
- Relationship between _endogenous_ cortisol and bone needs more research [[3]](https://doi.org/10.3389/fendo.2018.00526)
- Estradiol levels cycles during premenopausal years, but this had no effect compared to mean cycling level. It is suggested that cortisol is higher in the follicular phase versus the luteal phase [[4]](https://doi.org/10.3389/fendo.2020.00311), but more research is needed regarding this topic as well. Preferrably day-to-day cortisol levels as estradiol levels vary a lot within follicular and lutal phase as well.
- If modeling a premenopausal state, other important hormones such as progesterone, luteinizing hormone (LH) and follicular-stimulating hormone (FSH) should be included.
- If modeling an oral contraceptive state, effects of different progestins and dosages should be taken into account as well as degree of suppression of ovarian hormones.
- Receptor activity (e.g. GR, different ERs) and enzyme activity (e.g. 11-beta-HSD1 and 11-beta-HSD2 and aromatase) would be nice to include in the model
- Since estradiol is synthesized from testosterone and estrone, this would be interesting to include.
- Ethinyl estradiol actions on hepatic protein production besides cortisol binding globulin (CBG). (Especially IGF-1 in relation to bone)
- Free vs bound hormones should be taken into account

## Sources
- [1]   Peterson, M. C., & Riggs, M. M. (2010). A physiologically based mathematical model of integrated calcium homeostasis and bone remodeling. In Bone (Vol. 46, Issue 1, pp. 49–63). Elsevier BV. https://doi.org/10.1016/j.bone.2009.08.053
- [2]   Karin, O., Raz, M., Tendler, A., Bar, A., Korem Kohanim, Y., Milo, T., & Alon, U. (2020). A new model for the HPA axis explains dysregulation of stress hormones on the timescale of weeks. In Molecular Systems Biology (Vol. 16, Issue 7). EMBO. https://doi.org/10.15252/msb.20209510
- [3]   Suarez-Bregua, P., Guerreiro, P. M., & Rotllant, J. (2018). Stress, Glucocorticoids and Bone: A Review From Mammals and Fish. In Frontiers in Endocrinology (Vol. 9). Frontiers Media SA. https://doi.org/10.3389/fendo.2018.00526
- [4]   Hamidovic, A., Karapetyan, K., Serdarevic, F., Choi, S. H., Eisenlohr-Moul, T., & Pinna, G. (2020). Higher Circulating Cortisol in the Follicular vs. Luteal Phase of the Menstrual Cycle: A Meta-Analysis. In Frontiers in Endocrinology (Vol. 11). Frontiers Media SA. https://doi.org/10.3389/fendo.2020.00311
