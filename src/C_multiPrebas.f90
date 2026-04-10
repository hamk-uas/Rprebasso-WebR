 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!subroutine bridging  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine multiPrebas(multiOut,nSites,nClimID,nLayers,maxYears,maxThin, &
    nYears,thinning,pCrobas,allSP,siteInfo, maxNlayers, &
    nThinning,fAPAR,initClearcut,fixBAinitClarcut,initCLcutRatio,ETSy,P0y, initVar,&
    weatherPRELES,DOY,pPRELES, soilC,pYasso,&
    pAWEN,weatherYasso,litterSize,soilCtot, &
    defaultThin,ClCut,energyCuts,clct_pars,dailyPRELES,yassoRun,multiEnergyWood, &
    tapioPars,thdPer,limPer,ftTapio,tTapio,GVout,thinInt, &
    flagFert,nYearsFert,mortMod,pECMmod,& !protect removed btw nYearsFert and mortModX, neither in prebas subroutine nor multiPrebas() R function
    layerPRELES,LUEtrees,LUEgv, siteInfoDist, outDist, prebasFlags, &
	latitude, TsumSBBs)


implicit none

integer, parameter :: nVar=54,npar=53!, nSp=3
integer, intent(in) :: nSites, maxYears,maxThin,nClimID,maxNlayers,allSP
integer, intent(in) :: nYears(nSites),nLayers(nSites) !protect removed; neither in prebas subroutine nor multiPrebas() R function

 integer :: i,climID,ij,iz,ijj,ki,n,jj,az
 real (kind=8), intent(in) :: weatherPRELES(nClimID,maxYears,365,5)
 integer, intent(in) :: DOY(365),layerPRELES !, ECMmod fvec
 real (kind=8), intent(in) :: pPRELES(30),pCrobas(npar,allSP),tapioPars(5,2,3,20),pECMmod(12)
 real (kind=8), intent(inout) :: tTapio(5,allSP,2,7), ftTapio(5,allSP,3,7),mortMod(2)
 real (kind=8), intent(inout) :: siteInfo(nSites,11),thdPer(nSites),limPer(nSites)
 real (kind=8), intent(in) :: thinning(nSites,maxThin,11),pAWEN(12,allSP)
 real (kind=8), intent(inout) :: dailyPRELES(nSites,(maxYears*365),3)
 real (kind=8), intent(inout) :: LUEtrees(allSP),LUEgv
 real (kind=8), intent(inout) :: initClearcut(nSites,5),fixBAinitClarcut(nSites),initCLcutRatio(nSites,maxNlayers)  !initial stand conditions after clear cut. (H,D,totBA,Hc,Ainit)
! real (kind=8), intent(in) :: pSp1(npar),pSp2(npar),pSp3(npar)!,par_common
 real (kind=8), intent(in) :: defaultThin(nSites),ClCut(nSites),yassoRun(nSites)
 real (kind=8), intent(in) :: clct_pars(nSites,allSP,3),energyCuts(nSites)  !!energCuts
 real (kind=8), intent(in) :: thinInt(nSites) !site specific parameter that determines the thinning intensity; 
          !from below (thinInt>1) or above (thinInt<1);thinInt=999. uses the default value from tapio rules


! logical :: disturbanceON !!!this could be site specific but to block dist. in some sites you can work on the inputs
real (kind=8), intent(inout) :: siteInfoDist(nSites,10), outDist(nSites,maxYears,10) !inputs(siteInfoDist) & outputs(outDist) of disturbance modules

 !!! fertilization parameters
 !integer, intent(inout) :: fertThin !!! flag for implementing fertilization at thinning. the number can be used to indicate the type of thinning for now only thinning 3 fvec
 integer, intent(inout) :: flagFert !!! flag that indicates if fertilization has already been applied along the rotation
 integer, intent(inout) :: nYearsFert !!number of years for which the fertilization is effective

!!!ground vegetation
 !integer, intent(in) :: gvRun      !!!ground vegetation fvec
 real (kind=8), intent(inout) :: GVout(nSites,maxYears,5) !fAPAR_gv,litGV,photoGV,wGV      !!!ground vegetation
! integer, intent(in) :: siteThinning(nSites)
 integer, intent(inout) :: nThinning(nSites)
 real (kind=8), intent(out) :: fAPAR(nSites,maxYears)
 real (kind=8), intent(inout) :: initVar(nSites,7,maxNlayers),P0y(nClimID,maxYears,2),ETSy(nClimID,maxYears)!,par_common
 real (kind=8), intent(inout) :: multiOut(nSites,maxYears,nVar,maxNlayers,2),latitude(nSites), TsumSBBs(nSites,4)
 real (kind=8), intent(inout) :: multiEnergyWood(nSites,maxYears,maxNlayers,2)!!energCuts
 real (kind=8), intent(inout) :: soilC(nSites,maxYears,5,3,maxNlayers),soilCtot(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 ! real (kind=8) :: soilC(nSites,maxYears,5,3,maxNlayers),soilCtot(nSites,maxYears) !dimensions = nyears,AWENH,treeOrgans(woody,fineWoody,Foliage),species
 real (kind=8), intent(in) :: pYasso(35), weatherYasso(nClimID,maxYears,3),litterSize(3,allSP) !litterSize dimensions: treeOrgans,species
 real (kind=8) :: totBA(nSites), relBA(nSites,maxNlayers),mortModX
 real (kind=8) :: ClCutX, HarvArea,defaultThinX,maxState(nSites),check(maxYears)
 integer :: maxYearSite = 300,yearX(nSites),Ainit,sitex,ops(1),species

 integer :: etmodel,CO2model, gvRun, fertThin, ECMmod, oldLayer !not direct inputs anymore, but in prebasFlags fvec !wdimpl pflags
 integer, intent(inout) :: prebasFlags(10)

!!! 'un-vectorise' flags, fvec
etmodel = prebasFlags(1)
gvRun = prebasFlags(2)
fertThin = prebasFlags(3)
oldLayer = prebasFlags(4)
ECMmod = prebasFlags(5)
CO2model = prebasFlags(7)
! if(prebasFlags(6)==0) disturbanceON = .FALSE.
! if(prebasFlags(6)==1) disturbanceON = .TRUE.

!outDist(:,:,:) = 99!prebasFlags(6)
!outDist(1,10) = siteInfoDist(1,1)
!!!!initialize run
! multiOut = 0.
! open(1,file="test1.txt")
! open(2,file="test2.txt")

yearX = 0.
multiEnergyWood = 0.
!soilC = soilCinOut
!soilCtot = soilCtotInOut

! ---- initBiomasses loop: copy non-contiguous slices into contiguous temporaries ----
do i = 1,nSites
 do ijj = 1,nLayers(i)
  species = int(initVar(i,1,ijj))
  block
    real(kind=8) :: initVar_ib(7), biomasses_ib(nVar)
    initVar_ib(:) = initVar(i,:,ijj)
    biomasses_ib(:) = multiOut(i,1,:,ijj,1)
    call initBiomasses(pCrobas(:,species),initVar_ib,siteInfo(i,3),biomasses_ib,nVar,npar)
    initVar(i,:,ijj) = initVar_ib(:)
    multiOut(i,1,:,ijj,1) = biomasses_ib(:)
  end block
 enddo
enddo

! ---- main per-site prebas loop ----
do i = 1,nSites
 prebasFlags(8) = int(multiOut(i,1,7,1,2))
 multiOut(i,1,7,1,2) = 0.

  climID = siteInfo(i,2)
  defaultThinX = defaultThin(i)
  ClCutX = ClCut(i)
  
  !!!##set mortality model for managed and unmanaged forests
  mortModX = mortMod(1) !!mortality model to be used in the managed forests
  if(ClCut(i) < 0.5 .and. defaultThin(i) < 0.5) mortModX = mortMod(2) !!mortality model to be used in the unmanaged forests

  ! Copy non-contiguous per-site sections into contiguous work buffers before
  ! passing them to callees with explicit-shape dummy arrays. Some compilers
  ! materialize temporaries for these calls, but others do not, which can make
  ! the callee read the wrong stride. Large work arrays are heap-backed here to
  ! avoid exhausting the limited WebAssembly stack.
  block
    ! Small fixed-size temporaries.
    real(kind=8) :: siteInfo_tmp(11)
    real(kind=8) :: initClearcut_tmp(5)
    real(kind=8) :: siteInfoDist_tmp(10)
    real(kind=8) :: TsumSBBs_tmp(4)

    ! Small 2-D temporaries with bounded extents.
    real(kind=8) :: clct_pars_tmp(allSP, 3)

    ! Large variable-size temporaries are heap-backed to keep the Wasm stack small.
    real(kind=8), allocatable :: initCLcutRatio_tmp(:)
    real(kind=8), allocatable :: fAPAR_tmp(:)
    real(kind=8), allocatable :: ETSy_tmp(:)
    real(kind=8), allocatable :: soilCtot_tmp(:)
    real(kind=8), allocatable :: initVar_tmp(:,:)
    real(kind=8), allocatable :: thinning_tmp(:,:)
    real(kind=8), allocatable :: output_tmp(:,:,:,:)
    real(kind=8), allocatable :: P0y_tmp(:,:)
    real(kind=8), allocatable :: weatherPRELES_tmp(:,:,:)
    real(kind=8), allocatable :: soilC_tmp(:,:,:,:)
    real(kind=8), allocatable :: weatherYasso_tmp(:,:)
    real(kind=8), allocatable :: dailyPRELES_tmp(:,:)
    real(kind=8), allocatable :: energyWood_tmp(:,:,:)
    real(kind=8), allocatable :: GVout_tmp(:,:)
    real(kind=8), allocatable :: outDist_tmp(:,:)

    allocate(initCLcutRatio_tmp(maxNlayers))
    allocate(fAPAR_tmp(maxYears))
    allocate(ETSy_tmp(maxYears))
    allocate(soilCtot_tmp(maxYears))
    allocate(initVar_tmp(7, maxNlayers))
    allocate(thinning_tmp(nThinning(i), 11))
    allocate(output_tmp(nYears(i), nVar, nLayers(i), 2))
    allocate(P0y_tmp(nYears(i), 2))
    allocate(weatherPRELES_tmp(nYears(i), 365, 5))
    allocate(soilC_tmp(nYears(i), 5, 3, nLayers(i)))
    allocate(weatherYasso_tmp(nYears(i), 3))
    allocate(dailyPRELES_tmp(nYears(i)*365, 3))
    allocate(energyWood_tmp(nYears(i), nLayers(i), 2))
    allocate(GVout_tmp(nYears(i), 5))
    allocate(outDist_tmp(nYears(i), 10))

    ! ---- copy in ----
    siteInfo_tmp(:) = siteInfo(i,:)
    initVar_tmp(:, 1:nLayers(i)) = initVar(i,:,1:nLayers(i))
    thinning_tmp(:,:) = thinning(i, 1:nThinning(i), :)
    output_tmp(:,:,:,:) = multiOut(i, 1:nYears(i), :, 1:nLayers(i), :)
    fAPAR_tmp(1:nYears(i)) = fAPAR(i, 1:nYears(i))
    initClearcut_tmp(:) = initClearcut(i,:)
    initCLcutRatio_tmp(1:nLayers(i)) = initCLcutRatio(i, 1:nLayers(i))
    ETSy_tmp(1:nYears(i)) = ETSy(climID, 1:nYears(i))
    P0y_tmp(:,:) = P0y(climID, 1:nYears(i), :)
    weatherPRELES_tmp(:,:,:) = weatherPRELES(climID, 1:nYears(i), :, :)
    soilC_tmp(:,:,:,:) = soilC(i, 1:nYears(i), :, :, 1:nLayers(i))
    weatherYasso_tmp(:,:) = weatherYasso(climID, 1:nYears(i), :)
    soilCtot_tmp(1:nYears(i)) = soilCtot(i, 1:nYears(i))
    clct_pars_tmp(:,:) = clct_pars(i,:,:)
    dailyPRELES_tmp(:,:) = dailyPRELES(i, 1:(nYears(i)*365), :)
    energyWood_tmp(:,:,:) = multiEnergyWood(i, 1:nYears(i), 1:nLayers(i), :)
    GVout_tmp(:,:) = GVout(i, 1:nYears(i), :)
    siteInfoDist_tmp(:) = siteInfoDist(i,:)
    outDist_tmp(:,:) = outDist(i, 1:nYears(i), :)
    TsumSBBs_tmp(:) = TsumSBBs(i,:)

    call prebas(nYears(i),nLayers(i),allSP,siteInfo_tmp,pCrobas,initVar_tmp,&
    thinning_tmp,output_tmp,nThinning(i),maxYearSite,fAPAR_tmp, &
    initClearcut_tmp,fixBAinitClarcut(i),initCLcutRatio_tmp,ETSy_tmp,&
    P0y_tmp,weatherPRELES_tmp,DOY,pPRELES, &
    soilC_tmp,pYasso,pAWEN,weatherYasso_tmp,&
    litterSize,soilCtot_tmp,defaultThinX,&
    ClCutX,energyCuts(i),clct_pars_tmp,dailyPRELES_tmp,yassoRun(i),&
    energyWood_tmp,tapioPars,thdPer(i),limPer(i),ftTapio,tTapio,&
    GVout_tmp,thinInt(i), &
    flagFert,nYearsFert,mortModX,pECMmod,layerPRELES,LUEtrees,LUEgv, &
    siteInfoDist_tmp, outDist_tmp, prebasFlags,latitude(i), TsumSBBs_tmp)

    ! ---- copy back (inout arguments) ----
    siteInfo(i,:) = siteInfo_tmp(:)
    initVar(i,:,1:nLayers(i)) = initVar_tmp(:, 1:nLayers(i))
    multiOut(i, 1:nYears(i), :, 1:nLayers(i), :) = output_tmp(:,:,:,:)
    fAPAR(i, 1:nYears(i)) = fAPAR_tmp(1:nYears(i))
    initClearcut(i,:) = initClearcut_tmp(:)
    initCLcutRatio(i, 1:nLayers(i)) = initCLcutRatio_tmp(1:nLayers(i))
    ETSy(climID, 1:nYears(i)) = ETSy_tmp(1:nYears(i))
    P0y(climID, 1:nYears(i), :) = P0y_tmp(:,:)
    soilC(i, 1:nYears(i), :, :, 1:nLayers(i)) = soilC_tmp(:,:,:,:)
    soilCtot(i, 1:nYears(i)) = soilCtot_tmp(1:nYears(i))
    dailyPRELES(i, 1:(nYears(i)*365), :) = dailyPRELES_tmp(:,:)
    multiEnergyWood(i, 1:nYears(i), 1:nLayers(i), :) = energyWood_tmp(:,:,:)
    GVout(i, 1:nYears(i), :) = GVout_tmp(:,:)
    siteInfoDist(i,:) = siteInfoDist_tmp(:)
    outDist(i, 1:nYears(i), :) = outDist_tmp(:,:)
    TsumSBBs(i,:) = TsumSBBs_tmp(:)
  end block
end do
 ! close(1)
 ! close(2)
end subroutine

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
!subroutine bridging  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
