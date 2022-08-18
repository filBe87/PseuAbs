### =========================================================================
### define wsl.samplePseuAbs
### =========================================================================
#' Sample pseudoabsences using various strategies
#'
#' Flexlible function to sample pseudoabsences with various strategies and and thin
#' presences and pseudoabsences with flexible distance constraints. This is the
#' core function of the PseuAbs package.
#'
#' @param n Positive integer; number of pseudoabsence points desired. Default is 10000.
#' @param env.stack RasterStack/RasterBrick with environmental layers for sampling and
#' extraction
#' @param type Character; desired sampling strategy. Options are 'geographic', 'density',
#' 'random', 'target.group', 'geo.strat', 'env.strat' and 'env.semi.strat'
#' (see details). Default is 'geographic'.
#' @param add.strat Fraction between 0 and 1; should strategy be complemented
#' by a fraction of environmental strata? (Does not apply when type env.strat or
#' env.semi.strat is chosen)
#' @param pres SpatialPoints object; location of presence points. Necessary for
#' 'geographic' and 'density' strategies, otherwise optional. Presence
#' points can be thinned according to adhere to flexible distance constraints.
#' @param taxon Character (optional); name of taxon of interest to keep track of in meta
#' information.
#' @param geodist_fact Positive floating point number to adjust spatial autocorrelation lengths: for
#' 'geographic' pseudoabsence point patterns, values below 1 increase
#' autocorrelation length; values above 1 decrease it; for 'density' sampling
#' it is the other way around.
#' @param geores_fact Positive integer; aggregation factor for template raster from which pseudoabsences
#' are sampled (for 'geographic', 'density', 'random', and 'geo.strat' strategies). Larger
#' values save computation time, but derease resolution of sampling points.
#' @param template_dir Character; directory where template raster should be saved in/loaded from. Default is
#' tempdir(). The template will be saved in/loaded from the directory depending on
#' whether a template has already been created by a previous function call.
#' @param geo_nrep Positive Integer; number of replicates of models fitted for the 'geographic' strategy.
#' More will create a smoother pattern but increase computation time. Default is seven.
#' @param target.group_dir Character; directory where xy files of traget group taxa are
#' stored. Must be supplied if sampling stragegy is 'target.group', must contain
#' a column names 'x' and 'y' with coordinates in the same projection as the RasterStack/RasterBrick
#' provided to the env.stack argument.
#' @param env.strat_path Character; directory where sample of environmental strata for
#' 'env.strat' or 'env.semi.strat' sampling should be saved in/loaded from.
#' If NA, nothing will be saved; if provided, environmental strata will be saved in/loaded from
#' directory depending on whether a file has been created by a previous function call.
#' @param rAll Boolean; should all raster data be read into memory for computation of environmental
#' strata? This is faster but you may run into memory issues for large rasters. Default is TRUE.
#' @param force_spat_thin Character; should minimum distance be enforced between points? Options
#' are 'no', 'presences', 'absences', and 'both'. By default thinning is defined for
#' pseudoabsences from 'geographic', 'density', 'random', and 'geo.strat' methods
#' with minimum distance according to the resoultion of the template raster (see argument geores_fact).
#' 'presences' takes the minimum distance criterion from the template raster over
#' to the 'presence' points; 'absences' takes the criterion over to 'env.strat',
#' 'env.semi.strat', and 'target.group'; 'both' does it for both.
#' @param limdist Positive float; the minimum distance accepted for spatial thinning. Units
#' should be km if the spatial data is projected, otherwise the units of the
#' coordinate reference system are used. If no value is supplied, the maximum
#' distance between two cell centres of the template raster will be taken (see above).
#' @param set_max_npres_to_nabs Boolean; should the maximum number of presences be
#' equal to the number of pseudoabsences defined. Default is TRUE.
#' @details 'geographic' samples pseudoabsences with a sampling probability.
#' inversely proportional to the geographic distance to presence observations.
#' density' samples pseudoabsences proportional to the density of presence
#' observations. 'random' samples pseudoabsences randomly with a sampling
#' probability proportional the area of the cells. 'target.group' samples
#' pseudoabsences from the presences of the taxa of the target group, attempting
#' to correct for sampling bias. It depends on a directory with taxa defined by the
#' user as target group. 'geo.strat' samples pseudoabsences geographically
#' stratified either on a plane, or on a sphere depending on the projection of
#' the supplied env.stack. 'env.strat' samples pseudoabsences environmentally
#' stratified. Points are sampled form all realized combinations of environmental
#' conditions occuring in the environmental stack that have a minimal occurrence
#' frequency. Environmental strata are calculated based on all
#' raster layers supplied. If a directory is supplied as 'env.strat_path', a
#' large sample of stratified points will be saved to speed up computations for
#' follow-up species. If environmental strata based on different predictors than
#' supplied are preferred 'env.strat_path' can be an .RData file from a
#' previous sampling of strata form different environmental predictors. 'env.semi.strat'
#' is similar to 'env.strat' but samples environmental strata proportional to the
#' logarithm of the area they cover. See Descombes et al. for more details.
#'
#' @return a S4 object of class 'wsl.pseudoabsences'. It contains the follwing slots:
#' meta, a list with meta information; pa, a vector with 1 encoding presences and 0 encoding
#' pseudoabsences; xy, two-column matrix with coordinates for the (thinned) presence points
#' and the sampled pseudoabsences. The projection is the same as the one of the env.stack supplied; env_vars,
#' a data.frame containing extractions of the env.stack at the points xy; and call, which is
#' the call made to the wsl.samplePseuAbs function.
#' @author Philipp Brun
#' @export
#' @examples
#'### =========================================================================
#'### Data preparation
#'### =========================================================================
#'
#'# Predictors
#'bio=raster::getData('worldclim',var='bio',lon=16, lat=48,res=.5)
#'bio=bio[[c(1,4,12)]]
#'
#'# install.packages("rgbif")
#'library(rgbif)
#'# extract species
#'spn='Boletus aestivalis'
#'xt=as.vector(extent(bio))
#'baest <- occ_search(scientificName=spn,
#'                   hasCoordinate=TRUE,
#'                   decimalLongitude=paste0(xt[1],",",xt[3]),
#'                   decimalLatitude=paste0(xt[2],",",xt[4]))
#'
#'pbaest=baest$data[,c('decimalLongitude','decimalLatitude')]
#'baest_spp=SpatialPoints(pbaest,proj4string = crs(bio))
#'
#'# extract target group
#'targr <- occ_search(familyKey = 8789,
#'                    hasCoordinate=TRUE,
#'                    limit = 10000,
#'                    decimalLongitude=paste0(xt[1],",",xt[3]),
#'                    decimalLatitude=paste0(xt[2],",",xt[4]))
#'
#'ptargr=as.matrix(targr$data[,c('decimalLongitude','decimalLatitude')])
#'colnames(ptargr)=c("x","y")
#'
#'# create temporary directory for target.group info
#'tdir=paste0(tempdir(),"/trgr")
#'dir.create(tdir)
#'write.table(ptargr,file=paste0(tdir,"/targetxy.txt"),row.names = F)
#'
#'# create temporary directory for template raster and env strata
#'strdir=paste0(tempdir(),"/str")
#'dir.create(strdir)
#'
#'# Note that for these should not be temporary files for a real analysis.
#'
#'### =========================================================================
#'### Sample pseudoabsences
#'### =========================================================================
#'
#'# Geograhpic method with 20% env strata
#'pseu.abs1=wsl.samplePseuAbs(type="geographic",
#'                           n=5000,
#'                           env.stack=bio,
#'                           pres=baest_spp,
#'                           add.strat=0.2,
#'                           template_dir=strdir,
#'                           env.strat_path=strdir,
#'                           geodist_fact=1,
#'                           geores_fact=3,
#'                           geo_nrep=7,
#'                           taxon=spn)
#'
#'plot(pseu.abs1)
#'
#'# Only geographic with longer autocorrelation length
#'pseu.abs2=wsl.samplePseuAbs(type="geographic",
#'                            n=5000,
#'                            env.stack=bio,
#'                            pres=baest_spp,
#'                            add.strat=0,
#'                            template_dir=strdir,
#'                            env.strat_path=strdir,
#'                            geodist_fact=.5,
#'                            geores_fact=3,
#'                            geo_nrep=7,
#'                            taxon=spn)
#'
#'plot(pseu.abs2)
#'
#'# Random and thin presences by default resolution (determined by
#'# template raster)
#'pseu.abs3=wsl.samplePseuAbs(type="random",
#'                            n=5000,
#'                            env.stack=bio,
#'                            template_dir=strdir,
#'                            pres=baest_spp,
#'                            geores_fact=3,
#'                            add.strat=0,
#'                            taxon=spn,
#'                            force_spat_thin="presences")
#'
#'plot(pseu.abs3)
#'
#'# Geo.start and thin presences by min 10 km distance
#'pseu.abs4=wsl.samplePseuAbs(type="geo.strat",
#'                            n=5000,
#'                            env.stack=bio,
#'                            template_dir=strdir,
#'                            pres=baest_spp,
#'                            geores_fact=3,
#'                            add.strat=0,
#'                            taxon=spn,
#'                            force_spat_thin="presences",
#'                            limdist=10)
#'
#'plot(pseu.abs4)
#'
#'# Target group with 20% env.strat & thining presences & absences
#'pseu.abs5=wsl.samplePseuAbs(type="target.group",
#'                          n=5000,
#'                          env.stack=bio,
#'                          template_dir=strdir,
#'                          target.group_dir=tdir,
#'                          env.strat_path=strdir,
#'                          geores_fact=3,
#'                          pres=baest_spp,
#'                          add.strat=0.2,
#'                          taxon=spn,
#'                          force_spat_thin="both")
#'
#'plot(pseu.abs5)
#'
#'# Environmental semi-stratified
# pseu.abs6=wsl.samplePseuAbs(n = 5000,
#                            env.stack=bio,
#                            type = "env.semi.strat",
#                            add.strat = 0,
#                            pres = baest_spp,
#                            taxon = spn,
#                            template_dir=strdir,
#                            env.strat_path=strdir)
#'
#'plot(pseu.abs6)
#'
#'# Environmental semi-stratified with a non-feasible threshold of
#'# 30km limdist. CAREFUL, this one takes a few minutes!
#'pseu.abs7=wsl.samplePseuAbs(n = 5000,
#'                            env.stack=bio,
#'                            type = "env.semi.strat",
#'                            add.strat = 0,
#'                            pres = baest_spp,
#'                            taxon = spn,
#'                            template_dir=strdir,
#'                            env.strat_path=strdir,
#'                            force_spat_thin="both",
#'                            limdist=30)
#'
#'plot(pseu_abs7)
#'
#'# Density dependent
#'pseu.abs8=wsl.samplePseuAbs(n = 5000,
#'                            env.stack=bio,
#'                            type = "density",
#'                            add.strat = 0,
#'                            pres = baest_spp,
#'                            taxon = spn,
#'                            geores_fact=3,
#'                            template_dir=strdir,
#'                            env.strat_path=strdir)
#'
#'plot(pseu.abs8)
#'
wsl.samplePseuAbs<-function(n=10000,
                            env.stack,
                            type="geographic",
                            add.strat=0,
                            pres=numeric(),
                            taxon=character(),
                            geodist_fact=1,
                            geores_fact=1,
                            template_dir=tempdir(),
                            geo_nrep=7,
                            target.group_dir=NA,
                            env.strat_path=NA,
                            rAll=TRUE,
                            force_spat_thin="no",
                            limdist=NA,
                            set_max_npres_to_nabs=TRUE){

  ### ------------------------
  ### Check input and prepare
  ### ------------------------

  if(add.strat<0 | add.strat>1){
    stop("add.strat represents the fraction of pseudoabsences that are
         sampled environmentally stratified and should be between 0 and 1!")
  }

  possibtype=c("geographic","random","target.group","geo.strat","density",
               "env.strat","env.semi.strat")
  if(length(type)!=1 | !(type%in%possibtype)){
    stop("Invalid specification of pseudoabsence sampling type!")
  }

  possibthin=c("no","presences","absences","both")
  if(length(type)!=1 | !(force_spat_thin%in%possibthin)){
    stop("Invalid specification of spatial thinning method!")
  }

  if(grepl("env",type)){
    add.strat=1
  }

  ### ------------------------
  ### generate wsl.pseudoabsences object and add meta info
  ### ------------------------

  out<-preva.meta(type="pseudoabsence")

  tpnam=type
  if(tpnam=="geographic"){
    tpnam=paste0(tpnam,"_w",geodist_fact)
  }
  if(add.strat>0 && !grepl("env",type)){
    tpnam=paste0(tpnam,"_x_",add.strat,"env_strata")
  }

  out@meta$type=tpnam
  out@meta$taxon=taxon
  out@meta$force_spat_thin=force_spat_thin
  out@meta$template_file=NA

  call=match.call()
  out@call<-call

  ### ------------------------
  ### Prepare template raster
  ### ------------------------

  if(type%in%c("geographic","random","geo.strat","density") | force_spat_thin%in%c("presences","both")){

    # no template directory is supplied, just
    # calculate from scratch
    if(is.na(template_dir)){

      rst=aggregate(env.stack[[1]],
                    fact=geores_fact,
                    fun=function(x,na.rm){
                      ifelse(all(is.na(x)),NA,1)
                    },na.rm=T)
      # if template directory is supplied load if file exists
      # otherwise writeRaster
    } else {
      ptrn=paste0("template",geores_fact,".tif")
      tmfl=list.files(template_dir,pattern=ptrn,full.names=TRUE)

      if(length(tmfl)>0){
        rst=raster(tmfl[1])
      } else {
        tmfl=paste0(template_dir,"/",ptrn)
        if(geores_fact==1){
          rst=!is.na(env.stack[[1]])
          writeRaster(rst,filename=tmfl,overwrite=TRUE)
        } else {
          rst=aggregate(env.stack[[1]],
                        fact=geores_fact,
                        fun=function(x,na.rm){
                          ifelse(all(is.na(x)),NA,1)
                        },na.rm=T,filename=tmfl,overwrite=TRUE)
        }
      }
    }
    crs(rst)=crs(env.stack)

    # Calculate (latitudinal) distance between cells for potential
    # downstream analyses if no minim is provided
    if(is.na(limdist)){
      proje=grepl("longlat",crs(rst))
      dpp=SpatialPoints(coordinates(rst)[c(1,1+dim(rst)[2]),],
                        proj4string = crs(rst))
      limdist=spDists(dpp,longlat = proje)[1,2]
    }


    # Write path to template file into meta information
    out@meta$template_file=paste0(template_dir,"/template",geores_fact,".tif")
  }

  ### ------------------------
  ### Extract and refine presences
  ### ------------------------

  if(length(pres)>0){

    xt_pres=extract(env.stack,pres)
    sna=apply(xt_pres,1,function(x){
      any(is.na(x))
    })

    if(length(which(sna))>0){
      pres=pres[-which(sna)]
      xt_pres=xt_pres[-which(sna),]
      cat(paste0(length(which(sna))," non-matching presences removed..\n"))
    }

    if(force_spat_thin%in%c("presences","both")){

      n_mx=ifelse(n=="auto",50000,n)

      pres=SpatialPointsDataFrame(pres,data=data.frame(ID=1:length(pres)))
      if(nrow(pres@coords)<3000){
        tpp=thin_them(pres,limdist)
      } else {
        if(set_max_npres_to_nabs & nrow(pres@coords)<7*n_mx & nrow(pres@coords) >= n_mx){
          tpp=upsample_strategic(pres,limdist,n_mx,warnig=FALSE)
        }else if(set_max_npres_to_nabs & nrow(pres@coords)>=7*n_mx){
          tpp=upsample_thin(pres,limdist,n_mx)
        }else{
          tpp=upsample_strategic(pres,limdist,nrow(pres@coords),warnig=FALSE)
        }
      }

      cat(paste0(length(pres)-length(tpp)," presences removed to obtain min distance of ",
                 round(limdist,digits=2),"..\n"))
      pres=tpp
      xt_pres=xt_pres[pres$ID,]

      if(n=="auto"){
        n=ifelse(length(pres)<500,5000,ifelse(length(pres)>=500 & length(pres)<5000,10*length(pres),50000))
      }

    } else {

      if(n=="auto"){
        n=ifelse(length(pres)<500,5000,ifelse(length(pres)>=500 & length(pres)<5000,10*length(pres),50000))
      }

      if(set_max_npres_to_nabs & length(pres)>n){
        pres=pres[sample(1:length(pres),n)]
      }
    }
  }

  ### ------------------------
  ### Do the geographic sampling
  ### ------------------------

  if(type=="geographic"){

    # Sample a regular grid of abence points with
    # n_presences x geodist_fact points
    abs=sampleRegular(rst,round(length(pres)*geodist_fact),sp=T)

    # create geo absences from geo_nrep times jittering regular samples
    geo.pts=list()
    for(i in 1:geo_nrep){

      pt.abs=abs@coords[,c("x","y")]
      pt.abs=apply(pt.abs,2,jitter,factor=3)

      model.idw <- geoIDW(p=as.data.frame(pres@coords),
                          a=as.data.frame(pt.abs))

      prd <- predict(rst, model.idw,mask=TRUE)

      # Sample cell centers proportional to interpolated presence
      # probability
      nonaval=which(!is.na(values(prd)))

      prb=round(values(prd)[nonaval]*10^4)

      smp=sample(1:length(prb),size=round(n*(1-add.strat+.1)/geo_nrep),prob=prb,replace=F)

      geo.pts[[i]]=coordinates(prd)[nonaval[smp],]
    }

    # combine
    df.pseu=SpatialPoints(do.call("rbind",geo.pts),
                          proj4string = crs(rst))

    # extract and subsample to match desired number
    sp_abs=as(df.pseu,"SpatialPointsDataFrame")
    sp_abs@data=as.data.frame(extract(env.stack,df.pseu))

    ### ------------------------
    ### Do the random sampling
    ### ------------------------
  } else if (type=="random"){

    nona=which(!is.na(values(rst)))

    rnd.pts=sample(nona,
                   size=n*(1-add.strat)*1.5,
                   prob=values(suppressWarnings(raster::area(rst)))[nona])
    crds=coordinates(rst)[rnd.pts,]
    sp_abs=SpatialPointsDataFrame(coords=crds,
                                  data=as.data.frame(extract(env.stack,crds)),
                                  proj4string =rst@crs)

    ### ------------------------
    ### Do the target.group sampling
    ### ------------------------
  } else if (type=="target.group"){

    if(is.na("target.group_dir")){
      stop("target.group_dir has to be supplied for sampling type target.group!")
    }

    fls=list.files(target.group_dir,full.names=T)
    ltarpt=lapply(fls,read.table,header=TRUE)
    dftarpt=do.call("rbind",ltarpt)
    sptarpt=SpatialPoints(dftarpt,
                          proj4string = crs(env.stack))

    if(force_spat_thin%in%c("absences","both")){

      if(nrow(dftarpt)>7*n*(1-add.strat)*1.1){

        crds=upsample_thin(sptarpt,limdist,n*(1-add.strat)*1.1)

      } else {

        crds=upsample_strategic(sptarpt,limdist,n*(1-add.strat)*1.1)
      }
    } else {
      if(nrow(dftarpt)>n*(1-add.strat)*1.1){
        crds=sptarpt[sample(1:nrow(dftarpt),n*(1-add.strat)*1.1,replace=FALSE),]

      } else {
        crds=sptarpt[sample(1:nrow(dftarpt),n*(1-add.strat)*1.1,replace=TRUE),]
        cat(paste("Less target.group points available than requested. Sampling with replacemnent..."))

      }

    }

    sp_abs=SpatialPointsDataFrame(coords=crds,
                                  data=as.data.frame(extract(env.stack,crds)),
                                  proj4string =rst@crs)

    ### ------------------------
    ### Do the geo.strat sampling
    ### ------------------------
  } else if (type=="geo.strat"){

    if(grepl("longlat",rst@crs)){

      # sample regularly on the surface of a sphere
      fglb=(extent(rst)@xmax-extent(rst)@xmin)*(extent(rst)@ymax-extent(rst)@ymin)/(360*180)
      N=round(1.1*n*(1-add.strat)*length(values(rst))/length(which(!is.na(values(rst))))/fglb)
      r=1
      Nc=0
      a=4*pi*r^2/N
      d=sqrt(a)
      Mx=round(pi/d)
      dx=pi/Mx
      dy=a/dx
      pts=matrix(NA,ncol=3,nrow=N*5)

      for(m in 0:(Mx-1)){
        xx=pi*(m+0.5)/Mx
        My=round(2*pi*sin(xx)/dy)
        for(nn in 0:(My-1)){
          yy=2*pi*nn/My
          pts[Nc+1,]=c(r*sin(xx)*cos(yy),r*sin(xx)*sin(yy),r*cos(xx))
          Nc=Nc+1
        }
      }
      pts=na.omit(pts)

      # Transform to longitude latitude
      lat=atan(sqrt(pts[,2]^2+pts[,1]^2)/pts[,3])
      lat=ifelse(lat<0,lat+max(lat,na.rm = TRUE),lat-max(lat,na.rm=TRUE))
      lon=atan(pts[,2]/pts[,1])
      regpt=na.omit(unique(cbind(lon,lat)))
      regpt[,1]=regpt[,1]/pi*360
      regpt[,2]=regpt[,2]/pi*180
      sppt=SpatialPoints(regpt,proj4string = crs(rst))

      # crude extraction
      xtr1=extract(rst,sppt)
      sppt=sppt[which(!is.na(xtr1)),]

    } else {
      Ntarg=round(1.1*n*ncell(rst)/length(which(!is.na(values(rst)))))
      sppt=sampleRegular(rst,Ntarg,sp=T)
      sppt=sppt[-which(is.na(sppt@data[,1])),]

    }

    # Preprare output
    sp_abs=SpatialPointsDataFrame(coords=sppt@coords,
                                  data=as.data.frame(extract(env.stack,sppt)),
                                  proj4string =rst@crs)

  } else if (type=="density"){

    ### ------------------------
    ### Prepare point pattern object
    ### ------------------------

    xrng=extent(rst)[1:2]
    yrng=extent(rst)[3:4]

    # Define Point Pattern object to calculate
    owi=owin(xrange=xrng,yrange=yrng)
    myppp=ppp(x=pres@coords[,1],y=pres@coords[,2],window = owi)

    ### ------------------------
    ### Generate 'im' object with density info
    ### ------------------------

    lo=dim(rst)[1:2]

    x=seq(xrng[1],xrng[2],length.out = lo[2])
    y=seq(yrng[1],yrng[2],length.out = lo[1])

    dens=density(myppp,xy=list(x=x,y=y),adjust=geodist_fact/10)

    ### ------------------------
    ### Draw locations proportional to point density
    ### ------------------------

    drst=raster(dens)
    drst=resample(drst,rst)
    drst=mask(drst,rst)

    vls=values(drst)

    # Replace NA's with zero probability
    if(any(is.na(vls)) || any(vls<0)){
      vls[which(is.na(vls) | vls<0)]=0
    }

    # Sample from density distributions
    pts=sample(1:length(vls),n*(1-add.strat)*1.1,prob=vls)

    sppt=SpatialPoints(coordinates(drst)[pts,])

    # Preprare output
    sp_abs=SpatialPointsDataFrame(coords=sppt@coords,
                                  data=as.data.frame(extract(env.stack,sppt)),
                                  proj4string =rst@crs)

  }

  if(exists("sp_abs")){
    # Subsample and remove NAs
    nosna=apply(sp_abs@data,1,function(x){
      all(!is.na(x))
    })

    sp_abs=sp_abs[which(nosna),]
    if(nrow(sp_abs)>round(n*(1-add.strat))){
      sp_abs=sp_abs[sample(1:nrow(sp_abs),round(n*(1-add.strat))),]
    }

  }

  ### ------------------------
  ### Do the env.strat sampling
  ### ------------------------

  if (grepl("env",type) | add.strat>0){

    # define env strat sampling strategy
    if(add.strat<1){
      tyyp="env.strat"
    } else {
      tyyp=type
    }

    # if no strata directory is supplied, just
    # calculate from scratch
    if(is.na(env.strat_path)){

      strpts = create_envstrat(env.stk=env.stack,
                               rAll=rAll,
                               save_it=FALSE,
                               strat_dir=NA,
                               type=tyyp)

      # if strat directory is supplied load if file exists
      # otherwise write file
      # if .RData file is supplied read directly
    } else if(grepl(".RData",env.strat_path)){
      load(tmfl)
    } else{

      tmfl=list.files(env.strat_path,pattern=tyyp,full.names=TRUE)
      lyn=paste(names(env.stack),collapse="_")
      tmfl=grep(lyn,tmfl,value=TRUE)

      if(length(tmfl)>0){

        load(tmfl)
        strpts=strpts[sample(1:nrow(strpts)),]

      } else {

        strpts = create_envstrat(env.stk=env.stack,
                                 rAll=rAll,
                                 save_it=TRUE,
                                 strat_dir=env.strat_path,
                                 type=tyyp)

      }
    }

    if(force_spat_thin%in%c("absences","both")){

      strfy=upsample_thin(strpts,limdist,n*add.strat*1.1)

    } else {
      if(type=="env.semi.strat"){
        strfy=strpts[sample(1:nrow(strpts),n*add.strat*1.1),]
      } else {
        strfy=stratify(strpts,tyyp,n*add.strat*1.1)
      }
    }

    sp_abse=SpatialPointsDataFrame(coords=coordinates(strfy),
                                   data=as.data.frame(extract(env.stack,strfy)),
                                   proj4string =strpts@proj4string)

    nosna=apply(sp_abse@data,1,function(x){
      all(!is.na(x))
    })

    sp_abse=sp_abse[which(nosna),]
    sp_abse=sp_abse[sample(1:nrow(sp_abse),round(min(n*add.strat,nrow(sp_abse)))),]

    if(grepl("env",type)){
      sp_abs=sp_abse
    } else {
      sp_abs=rbind(sp_abse,sp_abs)
    }

  }

  ### ------------------------
  ### Prepare output
  ### ------------------------

  if(length(pres)>0){
    ptout=rbind(as(pres,"SpatialPoints"),as(sp_abs,"SpatialPoints"))

    out@pa=c(rep(1,length(pres)),
             rep(0,nrow(sp_abs)))

    out@env_vars=as.data.frame(rbind(xt_pres,sp_abs@data))
    out@xy=rbind(coordinates(pres),coordinates(sp_abs))
  } else {
    ptout=as(sp_abs,"SpatialPoints")
    out@pa=rep(0,nrow(sp_abs))
    out@env_vars=sp_abs@data
    out@xy=coordinates(sp_abs)
  }

  # return
  return(out)

}
