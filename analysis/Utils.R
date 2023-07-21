## Wrapper functions for repeated data.table procedures in the data management functions
## Set to global environment updates to avoid confusion between execution and global environments


#' Helper function to write data.table updates by reference
#'
#' Required to allow a data.table name passed as a string name
#' Otherwise unable to pass data.table as reference, alias, or character to function
#'
#' Function as suggested by Matt Dowle on SO to enable data.table to automatically optimise, as using get/mget or aliases prevents optimisations
#'
#' \describe{
#' \item{paste0}{Takes comma separated list of strings  and variable names containing strings NB does not need to be defined as a list or c}
#' \item{parse}{Converts the list to R call}
#' \item{eval}{Evaluates the R call within the parent environment}
#' }
#'
#' @return Nothing in itself, it constructs and evaluates a command in the parent environment
#'
EVAL <- function(...)eval(parse(text = paste0(...)), envir = .GlobalEnv)

#' Maximum by in place roll
#'
#' 7 times faster than using max and by
#' Superceded by max.grp.col_ which is faster and works across columns
#'
#' @param dt Data table name to update as char cannot be dt or then it does copy
#' @param aggregate.var Name char to apply min to
#' @param max.var Name char for group minimum
#' @param group Char vector of grouping variable/s
#'
#' @return  assign.dt assigned into parent frame
#'
max.roll_ <- function(dt,aggregate.var,max.var,group)  EVAL(dt,"[,",max.var, ":=",dt,"[",dt,"[,.I[which.max(",aggregate.var,")], keyby = ",group,"]$V1,c(",group,",'",aggregate.var,"')][",dt,"[,",group,"],",aggregate.var,", on = ",group,"]]")


#' Minimum by in place roll
#'
#' 7 times faster than using min and by
#'
#' @param dt Data table name to update as char cannot be dt or then it does copy
#' @param aggregate.var Name char to apply min to
#' @param min.var Name char for group minimum
#' @param group Char vector of grouping variable/s
#'
#' @return  assign.dt assigned into parent frame
#'
min.roll_ <- function(dt,aggregate.var,min.var,group) EVAL(dt,"[,",min.var, ":=",dt,"[",dt,"[,.I[which.min(",aggregate.var,")], keyby = ",group,"]$V1,c(",group,",'",aggregate.var,"')][",dt,"[,",group,"],",aggregate.var,", on = ",group,"]]")#, envir = .GlobalEnv)


#' Maximum by in place sort and roll across columns
#'
#' Note missing values are dropped by design before sort
#'
#' @param dt Data table name to update as char cannot be dt or then it does copy
#' @param aggregate.cols Name char to apply max to
#' @param max.var.name Name char for group max
#' @param id.vars Char vector of grouping variable/s
#'
#' @return  assign.dt assigned into parent frame
#'
max.grp.col_ <- function(dt, max.var.name, aggregate.cols, id.vars) { 
  EVAL(dt,'[,',max.var.name,' := data.table::melt(',
       dt,'[,c("',paste(c(id.vars,aggregate.cols),collapse = '","'),'"), with = F], na.rm = T, id.vars=c("'
       ,paste(c(id.vars),collapse = '","'),'"))[order(',paste(c(id.vars),collapse = ','),',-value),.SD[1], keyby = .(',paste(c(id.vars),collapse = ','),')][J(',dt,'[,.(',paste(c(id.vars),collapse = ','),')]),.(value)]]') }

#' Minimum by in place sort and roll across columns
#' 
#' Note missing values are dropped by design before sort
#'
#' @param dt Data table name to update as char cannot be dt or then it does copy
#' @param aggregate.cols Name char to apply min to
#' @param min.var.name Name char for group min
#' @param id.vars Char vector of grouping variable/s
#'
#' @return  assign.dt assigned into parent frame
#'
min.grp.col_ <- function(dt, min.var.name, aggregate.cols, id.vars) { 
  EVAL(dt,'[,',min.var.name,' := data.table::melt(',dt,'[,c("',paste(c(id.vars,aggregate.cols),collapse = '","'),'"), with = F], na.rm = T, id.vars=c("'
       ,paste(c(id.vars),collapse = '","'),'"))[order(',paste(c(id.vars),collapse = ','),',value),.SD[1], keyby = .(',paste(c(id.vars),collapse = ','),')][J(',dt,'[,.(',paste(c(id.vars),collapse = ','),')]),.(value)]]') }


#' Merge in place
#'
#' Expands patient datetimes to include all those in the merged data, so that a rolling left join can then be used
#'
#' @param dt Data table name to update as char cannot be dt or then it does copy
#' @param start.DTTM Date time variable name char in dt for start of row follow up time
#' @param end.DTTM Date time variable name char in dt for end of row follow up time
#' @param merge.dt Data table name char for merging in new data. This does have an additional variable start.DTTM added by functions
#' @param merge.start Date time variable char within merge.dt to align with start dates in dt
#' @param ID ID variable char defining the groups within which to calculate date sequences min / max date times.
#' @param assign.dt Name char for merged dataset to be assigned to. Can be the same as dt
#'
#' @return  assign.dt assigned into parent frame
#'
join.inplace_  <- function(dt,start.DTTM, end.DTTM, merge.dt, merge.start, ID, assign.dt) {
	# Sort by patient and date
	EVAL("data.table::setkey(",dt,",",ID,",",start.DTTM,")")
	EVAL(merge.dt,"[,", start.DTTM,":= ",merge.start,"]")
	# Sort mergeed by patient and date
	EVAL("data.table::setkey(",merge.dt,",",ID,",",merge.start,")")
	## Expand dataset rows to include all patient datetimes from merge
	# Done separately to ensure the minimum number of variables is used to save memory
	assign(x = assign.dt,
  			 value = EVAL(dt,"[data.table::SJ(unique(data.table::rbindlist(list(",dt,"[,c('",ID,"','",start.DTTM,"')],",
																																									merge.dt,"[,c('",ID,"','",start.DTTM,"')]",
				 						 "),use.names = T))),roll = Inf, rollends = c(T,T), mult = 'all']"),
  			 envir = .GlobalEnv)
	# Save latest end time as censor date per ID
	EVAL(assign.dt,"[,temp :=", end.DTTM,"[.N], by =",ID,"]")
	EVAL(assign.dt,"[,temp2 :=", start.DTTM,"[.N], by =",ID,"]")
	EVAL(assign.dt,"[temp2 > temp,temp :=temp2]")
	# Realign end date times
	EVAL(assign.dt,"[,", end.DTTM, ":= data.table::shift(", start.DTTM,", n = 1L, type = 'lead') - 1, by = ",ID,"]")
	# Final date should be the censor date and is any NA from previous row shift
	EVAL(assign.dt,"[is.na(",end.DTTM,"), ", end.DTTM,":= temp]")
	EVAL(assign.dt,"[, c('temp','temp2') := NULL]")
	# Finally perform the rolling join in place
	assign(x = assign.dt,
				 value = EVAL(assign.dt,"[, (names(",merge.dt,")[!(names(",merge.dt,") %in% c('",ID,"','",start.DTTM,"'))]) :=",
					 	merge.dt,"[data.table::SJ(",assign.dt,"[,c('",ID,"','",start.DTTM,"')]),, mult = 'all'][,c('",ID,"','",start.DTTM,"') := NULL]]"), envir = .GlobalEnv)
}



#' Extend sequential date range within group
#'
#' Expand range of dates by patient to match new data to allow left rolling joins
#' Existing unique date times in dt are kept, and new dates added
#' Range calculated by day from earliest date in range.dt
#'
#' @param dt Data table name to update as char cannot be dt or then it does copy
#' @param start.DTTM Date time variable name char in dt for start of row follow up time
#' @param end.DTTM Date time variable name char in dt for end of row follow up time
#' @param range.dt Data table name char where new date range to be derived from. NB this can be the same as dt for a self join.
#' @param range.start Date time variable char within range.dt to start sequence of dates. NB this can be the same as start.DTTM if using dt.
#' @param range.end Date time variable char within range.dt to end sequence of dates. NB this can be the same as start.DTTM if using dt.
#' @param ID ID variable char defining the groups within which to calculate date sequences min / max date times.
#'
#' @return  dt updated in parent frame
#'
dates.expand_  <- function(dt,start.DTTM, end.DTTM, range.dt, range.start, range.end, ID) {
	# Sort by patient and date
	EVAL("data.table::setkey(",dt,",",ID,",",start.DTTM,")")
	# Create combined unique dt of patients and date times across both datasets (can't do the following original anymore that repeated by day, as exceeds standard memory)
#	EVAL("indx <- unique(data.table::rbindlist(list(", range.dt,"[,.(", start.DTTM," = seq(from = min(",range.start,", na.rm = TRUE),
#                                                                                        to = max(",range.end,", na.rm = TRUE), by = 'day')),
#                                                                                        by = ",ID,"],",
	assign(x = 'indx',
				 value = EVAL("unique(data.table::rbindlist(list(",
			 range.dt,"[,min(",range.start,", na.rm = T), by = ",ID,"][,	data.table::setnames(.SD,2,'",start.DTTM,"')],",
			 dt,"[,.(",ID,",",start.DTTM,")],",
			 range.dt,"[,max(",range.end,", na.rm = T) - 1, by = ",ID,"][,	data.table::setnames(.SD,2,'",start.DTTM,"')])))"),
			 envir = .GlobalEnv)
	# Sort unique list of patients and date times
	EVAL("data.table::setkey(indx,",ID,",",start.DTTM,")")
	# Expand main dataset to all patient/ date time combinations using rolling join
	assign(x = dt,
				 value = EVAL(dt,"[indx,roll = T ,rollends = c(T,T)]"),
				 envir = .GlobalEnv)
	EVAL("rm(indx)")
	#rm(indx)
	# Sort by patient and date
	EVAL("data.table::setkey(",dt,",",ID,",",start.DTTM,")")
	# Realign end times for each row to additional start times
	EVAL(dt,"[,end := data.table::shift(",start.DTTM,", n = 1L, type = 'lead'), by = ",ID,"]")
	EVAL(dt,"[",ID," == data.table::shift(",ID,", n = 1L, type = 'lead'),",end.DTTM, ":= end - 1]")
	EVAL(dt,"[",ID," != data.table::shift(",ID,", n = 1L, type = 'lead'),",end.DTTM, ":= ",start.DTTM,"+ 1]")
	# Remove temporary variables
	EVAL(dt,"[,end := NULL]")
	# Remove patients in indx not in main dt (don't have an end.DTTM from main dt)
 EVAL(dt,"[!is.na(",end.DTTM,"),]")
	# Ensure data table is sorted
	EVAL("data.table::setkey(",dt,",",ID,",",start.DTTM,")")
#	assign(dt, value = EVAL(dt),
		#		 envir = .GlobalEnv)
}




#' Align start dates after a merge/roll
#'
#' Duplicate first row that a merged start date falls within, then align start and finish times for that row and subsequent set of matches within that row's time period
#'
#' @param dt Data table name to update as char -cannot be dt as then it does copy
#' @param start.DTTM Date time variable name char in dt for start of row follow up time
#' @param end.DTTM Date time variable name char in dt for end of row follow up time
#' @param merged.DTTM Date time variable name char in dt for that are the unaligned start merged date times
#' @param ID ID variable char defining the groups within which to align date times.
#'
#' @return  dt updated in parent frame
#'
#' @import data.table
#'
dates.expand.start.align_ <-  function(dt, start.DTTM, end.DTTM, ID, merged.DTTM) {
	# Sort by patient and date
	EVAL("data.table::setkey(",dt,",",ID,",",start.DTTM,",",merged.DTTM,")")
	# Duplicate first row in merged group where merged start date falls between start and end
	assign(x = dt, value = EVAL("data.table::rbindlist(list(",dt,",",dt,"[order(",ID,",",start.DTTM,",",merged.DTTM,")][",
								 merged.DTTM," > " ,start.DTTM, " & " , merged.DTTM , " < " ,end.DTTM, " & ",
								 start.DTTM , "!= data.table::shift(", start.DTTM ,", n = 1L, type = 'lag'),]))"), envir = .GlobalEnv)
	# Sort by patient and date
	EVAL("data.table::setkey(",dt,",",ID,",",start.DTTM,",",merged.DTTM,")")
	# Identify the new duplicated rows
	EVAL(dt,"[,GRP := .GRP, by = .(",ID,",",start.DTTM,",",merged.DTTM,")]")
	# Move 2nd row start date in group to merged start date
	EVAL(dt,"[", start.DTTM,"== data.table::shift(",start.DTTM,", n = 1L, type = 'lag'),",start.DTTM, ":=", merged.DTTM, "]")
	# Move 1st row end date in group to before merged start date
	EVAL(dt,"[GRP == data.table::shift(GRP, n = 1L, type = 'lead'),",end.DTTM, ":=", merged.DTTM, " - 0]")
	# For subsequent rows in matched sequence, identify time the subsequent row starts
	# (Done for all rows without filter as will be valid for all rows)
	EVAL(dt,"[,end := data.table::shift(",start.DTTM,", n = 1L, type = 'lead')]")
	# For subsequent rows in matched sequence set end to the time subsequent row ends
	# (Done for all rows without filter as will be valid for all rows)
	EVAL(dt,"[",ID," == data.table::shift(",ID,", n = 1L, type = 'lead'),",end.DTTM, ":= end - 0]")
	# Remove merged dates from outside matched sequence set
	EVAL(dt,"[GRP == data.table::shift(GRP, n = 1L, type = 'lead') & ",
			 merged.DTTM, "== data.table::shift(",merged.DTTM ,", n = 1L, type = 'lead'), ",merged.DTTM, ":= NA]")
	# Remove temporary variables
	EVAL(dt,"[,c('GRP','end') := NULL]")
	# Ensure data table is sorted
	EVAL("data.table::setkey(",dt,",",ID,",",start.DTTM,",",merged.DTTM,")")
#	assign(dt, value = EVAL(dt),
	#			 envir = .GlobalEnv)
}




#' Align end dates after a merge/roll
#'
#' Duplicate last row that a merged end date falls within, then align start and finish times for that matched row only
#'
#' @param dt Data table name to update as char - cannot be dt or then it does copy
#' @param start.DTTM Date time variable name char in dt for start of row follow up time
#' @param end.DTTM Date time variable name char in dt for end of row follow up time
#' @param merged.DTTM Date time variable name char in dt for that are the unaligned merged end date times
#' @param ID ID variable char defining the groups within which to align date times.
#'
#' @return  dt updated in parent frame
#'
#' @import data.table
#'
dates.expand.end.align_ <-  function(dt, start.DTTM, end.DTTM, ID, merged.DTTM) {
	# Sort by patient and date
	EVAL("data.table::setkey(",dt,",",ID,",",start.DTTM,",",merged.DTTM,")")
	# Duplicate last row in merged group where merged end date falls between start and end
	assign(dt,EVAL("data.table::rbindlist(list(",dt,",",dt,"[order(",ID,",",start.DTTM,",",merged.DTTM,")][",
								 merged.DTTM," > " ,start.DTTM, " & " , merged.DTTM , " < " ,end.DTTM, " & ",
								 end.DTTM , "!= data.table::shift(", end.DTTM ,",n = 1L, type = 'lead'),]))"),
				 , envir = .GlobalEnv)
	# Sort by patient and date
	EVAL("data.table::setkey(",dt,",",ID,",",start.DTTM,",",merged.DTTM,")")
	# Identify the new duplicated rows
	EVAL(dt,"[,GRP := .GRP, by = .(",ID,",",start.DTTM,",",merged.DTTM,")]")
	# Move 2nd row start date in group to after merged end date
	EVAL(dt,"[", start.DTTM,"==data.table::shift(",start.DTTM,", n = 1L, type = 'lag'),",start.DTTM, ":=", merged.DTTM ,"+ 0]")
	# Move 1st row end date in group to merged end date
	EVAL(dt,"[GRP == data.table::shift(GRP, n = 1L, type = 'lead'),",end.DTTM, ":=", merged.DTTM, "]")
	# For subsequent rows in matched sequence, no further updates are done, as should have been aligned first with dates.expand.start()
	# Remove merged dates from outside matched sequence set
	EVAL(dt,"[GRP == data.table::shift(GRP, n = 1L, type = 'lag') & ", merged.DTTM, "== data.table::shift(",merged.DTTM ,", n = 1L, type = 'lag'), ",merged.DTTM, ":= NA]")
	# Remove temporary variables
	EVAL(dt,"[,c('GRP') := NULL]")
	# Ensure data table is sorted
	EVAL("data.table::setkey(",dt,",",ID,",",start.DTTM,",",merged.DTTM,")")
#	assign(dt, value = EVAL(dt),
	#			 envir = .GlobalEnv)
}




#' Last observation carried forward
#'
#' Allows both numeric and character variables to be carried forward
#'
#' @param dt Data table name to update as char - cannot be dt or then it does copy
#' @param start.DTTM Date time variable name char in dt for start of row follow up time
#' @param ID ID variable char defining the groups within which to align date times.
#' @param group string list (surrounded by quotes) defining the change of group within which to limit locf
#' @param var.cols string list (surrounded by quotes) defining the variables to be carried forward
#'
#' @return  dt updated in parent frame
#'
#' @import data.table
#'
locf.roll_ <- function(dt, ID, start.DTTM, group, var.cols) {
	EVAL("if(sum(is.na(",dt,"[,",group,"]))>0) stop('missing in group')")
	EVAL(dt,"[, GRP := .GRP, by = ",group,"]")
	EVAL("data.table::setorder(",dt,",",ID,",",start.DTTM,",GRP)")
	EVAL("new_grp <- ",dt,"[order(",ID,",",start.DTTM,",GRP),
																c(TRUE, GRP[-1] != GRP[-.N])]")
	EVAL("new_grp[is.na(new_grp)] <- FALSE")  # If could be sure this wasn't needed then could condense function into one line?
	EVAL(dt,"[, (",var.cols,") :=  lapply(.SD,function(x) x[cummax((!is.na(x) | new_grp) * .I)]), .SDcols = ",var.cols,"]")
	EVAL("rm(new_grp)")
#	assign(dt, value = EVAL(dt),
#			 envir = .GlobalEnv)
}

nocb.roll_ <- function(dt, ID, start.DTTM, group, var.cols) {
  EVAL("if(sum(is.na(",dt,"[,",group,"]))>0) stop('missing in group')")
  EVAL(dt,"[, GRP := .GRP, by = ",group,"]")
  EVAL("data.table::setorder(",dt,",",ID,",-",start.DTTM,",GRP)")
  EVAL("new_grp <- ",dt,"[order(",ID,",-as.numeric(",start.DTTM,"),GRP),
																c(TRUE, GRP[-1] != GRP[-.N])]")
  EVAL("new_grp[is.na(new_grp)] <- FALSE")  # If could be sure this wasn't needed then could condense function into one line?
  EVAL(dt,"[, (",var.cols,") :=  lapply(.SD,function(x) x[cummax((!is.na(x) | new_grp) * .I)]), .SDcols = ",var.cols,"]")
  EVAL("rm(new_grp)")
  #	assign(dt, value = EVAL(dt),
  #			 envir = .GlobalEnv)
}

#' Delete rows in place for memory efficiency
#'  Taken from here https://github.com/Rdatatable/data.table/issues/635
#'
#' @import data.table
#'
del_rows <- function(X,delete) {

	keep <- !delete
	name_of_X <- deparse(substitute(X))
	X_names <- copy(names(X))
	X_new <- X[keep,X_names[1L],with=F]
	set(X,i=NULL,j=1L,value=NULL)

	for(j in seq_len(ncol(X))) {

		set(X_new,i=NULL,j=X_names[1L+j],value=X[[1L]][keep] )
		set(X,i=NULL,j=1L,value=NULL)

	}
 	assign(x = name_of_X,value=X_new, , envir = .GlobalEnv)
}


#' bootstrap p values for mixed effects linear model
#'
#' @import data.table
#' @import lme4
#'
lmer.fixef.glrt.bs <- function(iter, reduced.formula, model.data, full.formula, plot = F) {
	N <- iter
	boot.test.stats <- rep(0,N)

	fm.reduced <- lme4::lmer(formula = reduced.formula, data = model.data, REML = F)
	fm.full <- lme4::lmer(formula = full.formula, data = model.data, REML = F)
	obs.test.stat <- -2*(logLik(fm.reduced) - logLik(fm.full))
	attributes(obs.test.stat) <- NULL

	#1-pchisq(obs.test.stat,3)

	for(i in 1:N){
		model.data[,new.y := unlist(simulate(fm.reduced))]
		fm.reduced.new <- lme4::lmer(formula = update(reduced.formula, new.y ~ .), REML = F, data = model.data)
		fm.full.new <- lme4::lmer(formula = update(full.formula, new.y ~ .), REML = F, data = model.data)
		boot.test.stats[i] <- -2*(logLik(fm.reduced.new) - logLik(fm.full.new))
		if (plot == T) {
			par(mfrow = c(1,2))
			plot(as.numeric(model.data$SBP),model.data$new.y,ylim=c(5,20),pch=16,xaxp=c(1,4,3))
			title(paste("iteration:",i))
			hist(boot.test.stats[1:i],xlim=c(1,40),main="GLRT")
			points(boot.test.stats[i],0,pch=4,col="red", cex=2)
			abline(v=obs.test.stat, lwd=2)
		}
		#	if(i<10){readline()}
	}
	# Plot test statistics and chis-sq approximation
	if (plot == T) {
		par(mfrow=c(1,1))
		hist(boot.test.stats,prob=T,xlim=c(0,40),ylim=c(0,0.25))
		curve(dchisq(x,df=3),from=0,to=40,add=T,lwd=2)
		points(obs.test.stat,0,pch=4,col="red", cex=2)
	}

	# Get the p-value - how many simulated test statstics are larger than the observed one?
	# How does the p-value compare with the anova test p-value in step 3?
	return(mean(boot.test.stats>obs.test.stat))
}


#' bootstrap p values for mixed effects non linear model
#'
#' @import data.table
#' @import nlme
#'
nlmer.fixef.glrt.bs <- function(iter,
																reduced.formula, reduced.fixed, reduced.random, reduced.start, reduced.groups, reduced.formula.boot,
																model.data,
																full.formula, full.fixed, full.random, full.start, full.groups, full.formula.boot,
																plot = F) {
	N <- iter
	boot.test.stats <- rep(0,N)

	fm.reduced <- nlme::nlme(model = formula(reduced.formula), fixed = formula(reduced.fixed), random = formula(r1~1),
													 start = reduced.start, groups = formula(reduced.groups), data = model.data)
	fm.full <- nlme::nlme(model = formula(full.formula), fixed = formula(full.fixed), random = formula(full.random),
												start = full.start, groups = formula(full.groups), data = model.data)
	obs.test.stat <- -2*(logLik(fm.reduced) - logLik(fm.full))
	attributes(obs.test.stat) <- NULL

	#1-pchisq(obs.test.stat,3)
	full.patid <- model.data[,.(PATIENT_EXTERNAL_ID = as.character(PATIENT_EXTERNAL_ID))]
	for(i in 1:N){
		newparams <- MASS::mvrnorm(nrow(model.data), mu = nlme::fixed.effects(nl.model), Sigma = vcov(nl.model))
		randblups <- nlme::random.effects(nl.model)
		re <- data.table::data.table(PATIENT_EXTERNAL_ID = as.character(rownames(randblups)), randblups)

		model.data$new.y <- unlist(rep(1,nrow(model.data)) * newparams[1] +
			(model.data$CurrentSmoker * newparams[2]) +
			(newparams[3] * exp(-(newparams[4] + (newparams[5]-newparams[4])*model.data$CurrentSmoker ) * model.data$hour_inpatient)) +
			(model.data$male) * newparams[6] +
			(model.data$age51) * newparams[7] +
			(model.data$age61) * newparams[8] +
			(model.data$age71) * newparams[9] +
			(model.data$age81) * newparams[10] +
			(model.data$mixed) * newparams[11] +
			(model.data$Asian) * newparams[12] +
			(model.data$Black) * newparams[13] +
			(model.data$Other) * newparams[14] +
			(model.data$NotRecorded) * newparams[15] +
				re[full.patid, on = c('PATIENT_EXTERNAL_ID'), mult = 'all'][,2])

		fm.reduced.new <- nlme::nlme(model = reduced.formula.boot,  fixed = formula(reduced.fixed), random = formula(r1~1),
																 start = reduced.start, groups = formula(reduced.groups), data = model.data)
		fm.full.new <- nlme::nlme(model = full.formula.boot,  fixed = formula(full.fixed), random = formula(full.random),
															start = full.start, groups = formula(full.groups), data = model.data)
		boot.test.stats[i] <- -2*(logLik(fm.reduced.new) - logLik(fm.full.new))
		if (plot == T) {
			par(mfrow = c(1,2))
			plot(as.numeric(model.data$SBP),model.data$new.y,ylim=c(5,20),pch=16,xaxp=c(1,4,3))
			title(paste("iteration:",i))
			hist(boot.test.stats[1:i],xlim=c(1,40),main="GLRT")
			points(boot.test.stats[i],0,pch=4,col="red", cex=2)
			abline(v=obs.test.stat, lwd=2)
		}
		#	if(i<10){readline()}
	}
	# Plot test statistics and chis-sq approximation
	if (plot == T) {
		par(mfrow=c(1,1))
		hist(boot.test.stats,prob=T,xlim=c(0,40),ylim=c(0,0.25))
		curve(dchisq(x,df=3),from=0,to=40,add=T,lwd=2)
		points(obs.test.stat,0,pch=4,col="red", cex=2)
	}

	# Get the p-value - how many simulated test statstics are larger than the observed one?
	# How does the p-value compare with the anova test p-value in step 3?
	return(mean(boot.test.stats>obs.test.stat))
}



#' Convert factors to dummy variables within a data.table
#' Note does not exclude a baseline group
#'
#' @import data.table
#'
factor.to.dummy <- function(xdt, xcol) {
	EVAL('newcolnames <- ',xdt,'[,levels(as.factor(',xcol,'))][-1]')
	EVAL(xdt,'[,(newcolnames) := lapply(newcolnames, function(factor.col) as.numeric(',xcol,' == factor.col))]')
}


########################### Other functions for report that aren't data.table wrappers

numbers2words <- function(x){
  ## Function by John Fox found here: 
  ## http://tolstoy.newcastle.edu.au/R/help/05/04/2715.html
  ## Tweaks by AJH to add commas and "and"
  #https://github.com/ateucher/useful_code/blob/master/R/numbers2words.r
  #https://gist.github.com/psychemedia/numbers2words.R
  
  helper <- function(x){
    
    digits <- rev(strsplit(as.character(x), "")[[1]])
    nDigits <- length(digits)
    if (nDigits == 1) as.vector(ones[digits])
    else if (nDigits == 2)
      if (x <= 19) as.vector(teens[digits[1]])
    else trim(paste(tens[digits[2]],
                    Recall(as.numeric(digits[1]))))
    else if (nDigits == 3) trim(paste(ones[digits[3]], "hundred and", 
                                      Recall(makeNumber(digits[2:1]))))
    else {
      nSuffix <- ((nDigits + 2) %/% 3) - 1
      if (nSuffix > length(suffixes)) stop(paste(x, "is too large!"))
      trim(paste(Recall(makeNumber(digits[
        nDigits:(3*nSuffix + 1)])),
        suffixes[nSuffix],"," ,
        Recall(makeNumber(digits[(3*nSuffix):1]))))
    }
  }
  trim <- function(text){
    #Tidy leading/trailing whitespace, space before comma
    text=gsub("^\ ", "", gsub("\ *$", "", gsub("\ ,",",",text)))
    #Clear any trailing " and"
    text=gsub(" and$","",text)
    #Clear any trailing comma
    gsub("\ *,$","",text)
  }  
  makeNumber <- function(...) as.numeric(paste(..., collapse=""))     
  #Disable scientific notation
  opts <- options(scipen=100) 
  on.exit(options(opts)) 
  ones <- c("", "one", "two", "three", "four", "five", "six", "seven",
            "eight", "nine") 
  names(ones) <- 0:9 
  teens <- c("ten", "eleven", "twelve", "thirteen", "fourteen", "fifteen",
             "sixteen", " seventeen", "eighteen", "nineteen")
  names(teens) <- 0:9 
  tens <- c("twenty", "thirty", "forty", "fifty", "sixty", "seventy", "eighty",
            "ninety") 
  names(tens) <- 2:9 
  x <- round(x)
  suffixes <- c("thousand", "million", "billion", "trillion")     
  if (length(x) > 1) return(trim(sapply(x, helper)))
  helper(x)
}

rnd <- function(x) floor(x/10)*10

max.category <- function(predi) {unique(dt.tv[!is.na(get(covariates[predi])),get(covariates[predi])], na.rm = T)[which.max(dt.tv[!is.na(get(covariates[predi])),.N, by = c(covariates[predi])][['N']])]}

median.iqr <- function(x, dig=0) paste0(round(median(x, na.rm = T))," (",paste(round(quantile(x,c(0.25,0.75),na.rm = T),digits = dig),collapse = ","),")")

n.perc <- function(x, dig=0) paste0(rnd(sum(x == T, na.rm = T))," (",round(mean(x,na.rm = T),digits = dig)*100,"%)")

date.90.day <- function(x,ref.dat) is.finite(x) & as.numeric(x) - as.numeric(ref.dat) <= 90  & as.numeric(x) - as.numeric(ref.dat) >=0
date.30.day <- function(x,ref.dat) is.finite(x) & as.numeric(x) - as.numeric(ref.dat) <= 30  & as.numeric(x) - as.numeric(ref.dat) >=0

safelog <- function(x) { x[x < 1e-200] <- 1e-200; log(x) }

cuminc.km <- function(x,niter)  { 
  safelog <- function(x) { x[x < 1e-200] <- 1e-200; log(x) }
  n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) & !is.na(get(x)),event]))[-1]
  data.table::setkey(dt.tv, patient_id, tstart,tstop)
  est <- boot.est <- Reduce(function(x, y) merge(x, y, by = c("strata",'time'),all = T,sort = T, incomparables = 0),
                            lapply(n.type.events, function(i) data.table::setDT(summary(survival::survfit(survival::Surv(start,end,event==i) ~ get(x), 
                                                                                                          data = dt.tv[(postop.covid.cohort)  &  !is.na(get(x))], 
                                                                                                          id = patient_id),
                                                                                        times = sort(unique(dt.tv[(postop.covid.cohort) &    event == i & !is.na(get(x)),end])), extend = T
                            )[c('strata','time','n.risk', 'n.event')])[,(paste0('haz',i)) := n.event/n.risk][,c(1,2,5)])
  )[,(paste0('haz',n.type.events)) := lapply(.SD, function(x) data.table::fifelse(is.na(x),0,x)), .SDcols = paste0('haz',n.type.events)][
    order(strata,time),
    .(time,cumsum(exp(cumsum(log(1-Reduce('+',.SD))))*haz1)),
    keyby = strata, 
    .SDcols = paste0('haz',n.type.events)][time == 30,
                                           round(tail(V2,1),digits = 4)*100, 
                                           keyby = strata] 
  
  #  boot.est<- foreach::foreach(iter = 1:niter, .combine = 'cbind',.multicombine = T,.export = 'dt.tv',.packages = 'data.table') %dopar% {
  # 
  #   for(iter in 1:niter) {
  #     samp.patient_id <- dt.tv[patient_id %in% sort(sample(unique(dt.tv[(postop.covid.cohort) & !is.na(get(x)),patient_id]), replace = T))]
  #     data.table::setkey(dt.tv, patient_id, tstart,tstop)
  #      boot.est <- cbind(boot.est,Reduce(function(x, y) merge(x, y, by = c("strata",'time'),all = T,sort = T, incomparables = 0),
  #                                        lapply(n.type.events,
  #                                     function(i) {
  #        event.times = sort(unique(dt.tv[(postop.covid.cohort) & !is.na(get(x)) & event == i,end]))
  #        if(samp.patient_id[(event==i),.N] > 0) {
  #           data.table::setDT(summary(survival::survfit(survival::Surv(start,end,event==i) ~ get(x), 
  #       data = samp.patient_id,
  #       id = patient_id),
  #       times = event.times,
  #       extend = T
  #       )[c('strata','time','n.risk', 'n.event')])[,(paste0('haz',i)) := n.event/n.risk][,c(1,2,5)]
  #        }
  #        else {
  #         data.table::data.table('strata' = rep(levels(est[[1]]), 
  #                                               each = length(event.times)),
  #                                'time' = rep(event.times, 
  #                                             times = length(levels(est[[1]]))))[,(paste0('haz',i)) := 0] 
  #        }
  #      }
  #     )
  # )[,(paste0('haz',n.type.events)) := lapply(.SD, function(x) data.table::fifelse(is.na(x),0,x)),
  #   .SDcols = paste0('haz',n.type.events)][
  #   order(strata,time),.(time,
  #                        cumsum(exp(cumsum(safelog(1 - Reduce('+',.SD))))*haz1)),
  #   keyby = strata, 
  #   .SDcols = paste0('haz',n.type.events)][order(strata,time),
  #                                          round(tail(V2,1),digits = 3)*100,
  #                                          keyby = strata][,2])
  # }
  return(cbind(est[,1],
               rnd(data.table::setDT(summary(survival::survfit(survival::Surv(start,end,event==1) ~ get(x), 
                                                               data = dt.tv[(postop.covid.cohort) & !is.na(get(x)),],
                                                               id = patient_id), times = 0)[c('n.risk')])),
               rnd(data.table::setDT(summary(survival::survfit(survival::Surv(start,end,event==1) ~ get(x), 
                                                               data = dt.tv[(postop.covid.cohort) & !is.na(get(x)),],
                                                               id = patient_id), times = 30)[c('n.event')])),
               100*round(1 - (data.table::setDT(summary(survival::survfit(survival::Surv(start,end,event==1) ~ get(x), 
                                                                          data = dt.tv[(postop.covid.cohort) & !is.na(get(x)),],
                                                                          id = patient_id), times = 30)[c('surv')])), digits = 4),
               est[,2]
  )
  )
}


cuminc.km.sub <- function(x,niter)  { 
  safelog <- function(x) { x[x < 1e-200] <- 1e-200; log(x) }
  n.type.events <- sort(unique(dt.tv[(postop.covid.cohort) & sub.op == T & !is.na(get(x)),event]))[-1]
  data.table::setkey(dt.tv, patient_id, tstart,tstop)
  est <- boot.est <- Reduce(function(x, y) merge(x, y, by = c("strata",'time'),all = T,sort = T, incomparables = 0),
                            lapply(n.type.events, function(i) data.table::setDT(summary(survival::survfit(survival::Surv(start,end,event==i) ~ get(x), 
                                                                                                          data = dt.tv[(postop.covid.cohort) & sub.op == T & !is.na(get(x))], 
                                                                                                          id = patient_id), censored = F,
                                                                                        times = sort(unique(dt.tv[(postop.covid.cohort) & sub.op == T &    event == i & !is.na(get(x)),end])), 
                                                                                        extend = T
                            )[c('strata','time','n.risk', 'n.event')])[,(paste0('haz',i)) := n.event/n.risk][,c(1,2,5)])
  )[,(paste0('haz',n.type.events)) := lapply(.SD, function(x) data.table::fifelse(is.na(x),0,x)), .SDcols = paste0('haz',n.type.events)][
    order(strata,time),
    .(time,cumsum(exp(cumsum(safelog(1-Reduce('+',.SD))))*haz1)),
    keyby = strata, 
    .SDcols = paste0('haz',n.type.events)][time == 30,
                                           round(tail(V2,1),digits = 4)*100, 
                                           keyby = strata] 
  return(cbind(est[,1],
               rnd(data.table::setDT(summary(survival::survfit(survival::Surv(start,end,event==1) ~ get(x), 
                                                               data = dt.tv[(postop.covid.cohort) & !is.na(get(x)) & sub.op == T & start >= 0,],
                                                               id = patient_id), 
                                             times = 30,
                                             extend = T)[c('n.risk')])),
               rnd(data.table::setDT(summary(survival::survfit(survival::Surv(start,end,event==1) ~ get(x), 
                                                               data = dt.tv[(postop.covid.cohort) & !is.na(get(x)) & sub.op == T & start >= 0,],
                                                               id = patient_id), 
                                             times = 30,
                                             extend = T)[c('n.event')])),
               100*round(1 - (data.table::setDT(summary(survival::survfit(survival::Surv(start,end,event==1) ~ get(x), 
                                                                          data = dt.tv[(postop.covid.cohort) & !is.na(get(x)) & sub.op == T,],
                                                                          id = patient_id), times = 30,
                                                         extend = T)[c('surv')])), digits = 4),
               est[,2]
  )
  )
}



cuminc.km.mort <- function(x,niter)  { 
  safelog <- function(x) { x[x < 1e-200] <- 1e-200; log(x) }
  data.table::setkey(dt.tv, patient_id, tstart,tstop)

  return(cbind(rep(x,length(levels(dt.tv[start >=0 & any.op == T,
                                                     as.factor(get(x))]))),
               rnd(data.table::setDT(summary(survival::survfit(survival::Surv(start,end,died==1) ~ get(x), 
                                                               data = dt.tv[start >=0 & any.op == T & !is.na(get(x)),],
                                                               id = patient_id), times = 0)[c('n.risk')])),
               rnd(data.table::setDT(summary(survival::survfit(survival::Surv(start,end,died==1) ~ get(x), 
                                                               data = dt.tv[start >=0 & any.op == T & !is.na(get(x)),],
                                                               id = patient_id), times = 30)[c('n.event')])),
               100*round(1 - (data.table::setDT(summary(survival::survfit(survival::Surv(start,end,died==1) ~ get(x), 
                                                                          data = dt.tv[start >=0 & any.op == T & !is.na(get(x)),],
                                                                          id = patient_id), times = 30)[c('surv')])), digits = 4)
  )
  )
}

cuminc.km.mort.sub <- function(x,niter)  { 
  safelog <- function(x) { x[x < 1e-200] <- 1e-200; log(x) }
  data.table::setkey(dt.tv, patient_id, tstart,tstop)
  
  return(cbind(rep(x,length(levels(dt.tv[start >=0 & sub.op == T,
                                         as.factor(get(x))]))),
               rnd(data.table::setDT(summary(survival::survfit(survival::Surv(start,end,died==1) ~ get(x), 
                                                               data = dt.tv[start >=0 & sub.op == T & !is.na(get(x)),],
                                                               id = patient_id), times = 0)[c('n.risk')])),
               rnd(data.table::setDT(summary(survival::survfit(survival::Surv(start,end,died==1) ~ get(x), 
                                                               data = dt.tv[start >=0 & sub.op == T & !is.na(get(x)),],
                                                               id = patient_id), times = 30)[c('n.event')])),
               100*round(1 - (data.table::setDT(summary(survival::survfit(survival::Surv(start,end,died==1) ~ get(x), 
                                                                          data = dt.tv[start >=0 & sub.op == T & !is.na(get(x)),],
                                                                          id = patient_id), times = 30)[c('surv')])), digits = 4)
  )
  )
}

cuminc.cox <- function(n.type.events = c(1,2),dt, model, newdata, day) {
  eval(parse(text = paste0('assign(x = "base.haz", value = lapply(n.type.events, function(i) { survival::basehaz(',model,'[[i]],centered = F)[] }),envir=environment())')))
  eval(parse(text = paste0('assign(x = "base.haz", value = lapply(n.type.events, function(i) { base.haz[[i]][base.haz[[i]]$time %in% sort(unique(',dt,'[get(all.vars(get(model)[[i]]$call)[[3]])== i ,end])),][] }),envir=environment())')))
  base.haz.comp <- lapply(n.type.events, function(i) { data.table::data.table('time' = base.haz[[i]]$time,
                                                                         'base.haz' = base.haz[[i]][,1] - 
                                                                         c(0,head(base.haz[[i]][,1],-1)))})
  if(length(n.type.events) > 1) {
  base.haz.merge <- Reduce(x = base.haz.comp,f = function(x,y) merge(x,y,by = 'time', no.dups = T, suffixes = c(".x",".y"), all = T))
  } else {
    base.haz.merge <- base.haz.comp[[1]]
  }
  base.haz.merge[is.na(base.haz.merge)] <- 0
  
  eval(parse(text = paste0('assign(x = "risk",value = lapply(n.type.events, function(i) {
        data.table::data.table(
          "risk" = predict(object = ',model,'[[i]], 
                         type = "risk",
                         newdata = ',newdata,')) }),envir=environment())')))
  
  days.in.data <- vapply(1:length(day), FUN.VALUE = 1, FUN = function(i)  {
    ifelse(day[i] %in% unlist(base.haz.merge[,1]) == T,day[i],
     unlist(base.haz.merge[,1][max(which(base.haz.merge[,1] <= day[i])),1]))
  })
  
    return(matrix(100*round(apply(exp(-apply(Reduce('+',lapply(n.type.events, function(i) {
      outer(unlist(base.haz.merge[,.SD,.SDcols = (i+1)]) ,unlist(risk[[i]]),'*')})),2,cumsum)) *
        outer(unlist(base.haz.merge[,2]),unlist(risk[[1]]),'*'),2,cumsum), digits = 4), ncol = nrow(get(newdata)))[which(unlist(base.haz.merge[,1]) %in% days.in.data),])

}

# 
# new.data.postop.covid <- data.table::data.table('start.readmit' = 0,
#                                                 'end.readmit' = 0,
#                                                 'event.readmit' = F,
#                                                 'Abdominal' = F,
#                                                 'Cardiac'=F,
#                                                 'Obstetrics'=F,
#                                                 'Orthopaedic'=F,
#                                                 'Thoracic'=F,
#                                                 'Vascular'=F,
#                                                 'postcovid' =F,
#                                                 'age.cat' = levels(dt.tv$age.cat)[1],
#                                                 'sex' = levels(dt.tv$sex)[1],
#                                                 'bmi.cat' = levels(dt.tv$bmi.cat)[1],
#                                                 'imd5' = levels(dt.tv$imd5)[1],
#                                                 'wave' = 'Wave_1',
#                                                 'vaccination.status.factor' = levels(dt.tv$vaccination.status.factor)[1],
#                                                 'region' = levels(dt.tv$region)[1],
#                                                 'Current.Cancer' = F,
#                                                 'Emergency' =   F,
#                                                 'Charl12' =  levels(dt.tv$Charl12)[1],
#                                                 'recentCOVID' =F,
#                                                 'previousCOVID' = F,
#                                                 'patient_id' = 1)
# 
# Reduce(function(x, y) merge(x, y, by = c("strata",'time'),all = T,sort = T, incomparables = 0),
#        lapply(n.type.events, function(i) data.table::setDT(summary(survival::survfit(post.op.readmit.model[[i]], 
#                                                                                      data = dt.tv[(postop.covid.cohort)], 
#                                                                                      id = patient_id),
#                                                                    times = sort(unique(dt.tv[(postop.covid.cohort) & event == i ,end])), 
#                                                                    newdata = new.data.postop.covid,extend = T
#        )[c('cumhaz','surv','time','n.risk', 'n.event')]))
# )[,(paste0('haz',n.type.events)) := lapply(.SD, function(x) data.table::fifelse(is.na(x),0,x)), .SDcols = paste0('haz',n.type.events)][
#   order(strata,time),
#   .(time,cumsum(exp(cumsum(log(1-Reduce('+',.SD))))*haz1)),
#   keyby = strata, 
#   .SDcols = paste0('haz',n.type.events)][time == 30,
#                                          round(tail(V2,1),digits = 3)*100, 
#                                          keyby = strata] 


glance.flexsurvreg <- function(x,...) {
	tibble::tibble(
		logLik = x$loglik,
		AIC = x$AIC,
		df.residual = x$pars
	)
}

tidy.flexsurvreg <- function(x,...) {
	tibble::tibble(
		term = names(x$coefficients),
		estimate = x$coefficients,
		std.error = x$res[,4]
	)
}

#' Function to categorise ICD into ICD cause of death categories
#' 
#' fcase evaluates sequentially down list, So catch all at end of each chapter 
#' only used if no sub categories alreay matched
#' 
#' @import data.table
#' 
#' @param ICD,code ICD 10 code without puncuation
#' 
#' @return Character description of ICD 10 category
#' 
ICD_COD_categories <- function(ICD.code) {
  fcase(substr(ICD.code,1,2) == "A0" , "Intestinal infectious diseases",
        substr(ICD.code,1,3) %in% c("A15","A16") , "Respiratory tuberculosis",
        substr(ICD.code,1,3) %in% c("A17","A19") , "Other tuberculosis",
        substr(ICD.code,1,3) %in% c("A39") , "Meningococcal infection",
        substr(ICD.code,1,3) %in% c("A40","A41") , "Sepsis",
        substr(ICD.code,1,2) %in% c("B1") & as.numeric(substr(ICD.code,3,3)) >=5 , "Viral hepatitis",
        substr(ICD.code,1,2) %in% c("B2") & as.numeric(substr(ICD.code,3,3)) <5 , "HIV",
        substr(ICD.code,1,3) %in% c("B90"), "Sequelae of tuberculosis",
        substr(ICD.code,1,1) %in% c("B") , "Other infection",
        substr(ICD.code,1,2) %in% c("C0") , "Malignant neoplasms of lip, oral cavity and pharynx",
        substr(ICD.code,1,2) %in% c("C1") & as.numeric(substr(ICD.code,3,3)) <5 , "Malignant neoplasms of lip, oral cavity and pharynx",
        substr(ICD.code,1,3) %in% c("C15"), "Malignant neoplasm of oesophagus",
        substr(ICD.code,1,3) %in% c("C16"), "Malignant neoplasm of stomach",
        substr(ICD.code,1,3) %in% c("C18"), "Malignant neoplasm of colon",
        substr(ICD.code,1,3) %in% c("C19"), "Malignant neoplasm of rectosigmoid junction, rectum and anus",
        substr(ICD.code,1,2) %in% c("C2") & as.numeric(substr(ICD.code,3,3)) <=1 , "Malignant neoplasm of rectosigmoid junction, rectum and anus",
        substr(ICD.code,1,3) %in% c("C22"), "Malignant neoplasm of liver and intrahepatic bile ducts",
        substr(ICD.code,1,3) %in% c("C23","C24"), "Malignant neoplasm of gallbladder and biliary tract",
        substr(ICD.code,1,3) %in% c("C25"), "Malignant neoplasm of pancreas",
        substr(ICD.code,1,3) %in% c("C32"), "Malignant neoplasm of larynx",
        substr(ICD.code,1,3) %in% c("C33","C34"), "Malignant neoplasm of trachea, bronchus and lung",
        substr(ICD.code,1,3) %in% c("C43"), "Malignant melanoma of skin",
        substr(ICD.code,1,3) %in% c("C44"), "Other malignant neoplasms of skin",
        substr(ICD.code,1,3) %in% c("C45"), "Mesothelioma",
        substr(ICD.code,1,3) %in% c("C46"), "Kaposi sarcoma",
        substr(ICD.code,1,3) %in% c("C50"), "Malignant neoplasm of breast",
        substr(ICD.code,1,3) %in% c("C53"), "Malignant neoplasm of cervix uteri",
        substr(ICD.code,1,3) %in% c("C54","C55"), "Malignant neoplasm of other and unspecified parts of uterus",
        substr(ICD.code,1,3) %in% c("C56"), "Malignant neoplasm of ovary",
        substr(ICD.code,1,3) %in% c("C61"), "Malignant neoplasm of prostate",
        substr(ICD.code,1,3) %in% c("C62"), "Malignant neoplasm of testis",
        substr(ICD.code,1,3) %in% c("C64"), "Malignant neoplasm of kidney, except renal pelvis",
        substr(ICD.code,1,3) %in% c("C67"), "Malignant neoplasm of bladder",
        substr(ICD.code,1,3) %in% c("C71"), "Malignant neoplasm of brain",
        substr(ICD.code,1,3) %in% c("C81"), "Hodgkin lymphoma",
        substr(ICD.code,1,3) %in% c("C82","C83","C84","C85"), "Non-Hodgkin lymphoma",
        substr(ICD.code,1,3) %in% c("C90"), "Multiple myeloma and malignant plasma cell neoplasms",
        substr(ICD.code,1,3) %in% c("C91","C92","C93","C94","C95"), "Leukaemia",
        substr(ICD.code,1,3) %in% c("C97"), "Malignant neoplasms of independent (primary) multiple sites",
        substr(ICD.code,1,1) %in% c("C") , "Other malignancy",
        substr(ICD.code,1,1) %in% c("D") & as.numeric(substr(ICD.code,2,2)) <5 , "In situ and benign neoplasms, and neoplasms of uncertain or unknown behaviour",
        substr(ICD.code,1,1) %in% c("D") & as.numeric(substr(ICD.code,2,2)) %in% c(5,6)  , "Anaemias",
        substr(ICD.code,1,1) %in% c("E") & substr(ICD.code,2,2) == "1"  , "Diabetes",
        substr(ICD.code,1,3) %in% c("F01","F03"), "Vascular and unspecified dementia",
        substr(ICD.code,1,2) %in% c("F1"), "Mental and behavioural disorders due to psychoactive substance use",
        substr(ICD.code,1,1) %in% c("F") , "Other psychiatric",
        substr(ICD.code,1,3) %in% c("G00","G03"), "Meningitis (excluding meningococcal)",
        substr(ICD.code,1,4) %in% c("G122"), "Motor neuron disease",
        substr(ICD.code,1,3) %in% c("G20"), "Parkinson disease",
        substr(ICD.code,1,3) %in% c("G30"), "Alzheimer disease",
        substr(ICD.code,1,3) %in% c("G35"), "Multiple sclerosis",
        substr(ICD.code,1,3) %in% c("G40"), "Epilepsy",
        substr(ICD.code,1,1) %in% c("G") , "Other neurology",
        substr(ICD.code,1,3) %in% c("I21","I22"), "Acute myocardial infarction",
        substr(ICD.code,1,2) %in% c("I2") & as.numeric(substr(ICD.code,3,3)) >5 , "Other heart diseases",
        substr(ICD.code,1,1) %in% c("I") & substr(ICD.code,2,2) %in% c("3","4","5")  , "Other heart diseases",
        substr(ICD.code,1,3) %in% c("I60","I61","I62"), "Intracranial haemorrhage",
        substr(ICD.code,1,3) %in% c("I63"), "Cerebral infarction",
        substr(ICD.code,1,3) %in% c("I64"), "Stroke, not specified as haemorrhage or infarction",
        substr(ICD.code,1,3) %in% c("I70"), "Atherosclerosis",
        substr(ICD.code,1,3) %in% c("I71"), "Aortic aneurysm and dissection",
        substr(ICD.code,1,1) %in% c("I") , "Other cardiovascular",
        substr(ICD.code,1,3) %in% c("J09"), "Influenza due to certain identified influenza virus",
        substr(ICD.code,1,3) %in% c("J10","J11"), "Influenza",
        substr(ICD.code,1,2) %in% c("J1") & as.numeric(substr(ICD.code,3,3)) >=2  & as.numeric(substr(ICD.code,3,3)) <=8, "Pneumonia",
        substr(ICD.code,1,2) %in% c("J4") & as.numeric(substr(ICD.code,3,3)) <5, "Bronchitis, emphysema and other chronic obstructive pulmonary disease",
        substr(ICD.code,1,3) %in% c("J45","J46") , "Asthma",
        substr(ICD.code,1,1) %in% c("J") , "Other respiratory",
        substr(ICD.code,1,3) %in% c("K25","K26","K27") , "Gastric and duodenal ulcer",
        substr(ICD.code,1,2) %in% c("K4") & as.numeric(substr(ICD.code,3,3)) <=6, "Hernia",
        substr(ICD.code,1,3) %in% c("K57"), "Diverticular disease of intestine",
        substr(ICD.code,1,2) %in% c("K7") & as.numeric(substr(ICD.code,3,3)) <=7, "Diseases of the liver",
        substr(ICD.code,1,1) %in% c("K") , "Other GI",
        substr(ICD.code,1,1) %in% c("L") , "Diseases of the skin and subcutaneous tissue",
        substr(ICD.code,1,3) %in% c("M05","M06","M08") , "Rheumatoid arthritis and juvenile arthritis",
        substr(ICD.code,1,3) %in% c("M80","M81") , "	Osteoporosis",
        substr(ICD.code,1,1) %in% c("N") & as.numeric(substr(ICD.code,2,3)) <=15, "Glomerular and renal tubulo-interstitial diseases",
        substr(ICD.code,1,3) %in% c("N17","N18","N19") , "Renal failure",
        substr(ICD.code,1,4) %in% c("N390") , "UTI",
        substr(ICD.code,1,3) %in% c("R54") , "Hyperplasia of prostate",
        substr(ICD.code,1,3) %in% c("R95") , "Senility",
        substr(ICD.code,1,3) %in% c("R99") , "Sudden infant death syndrome",
        substr(ICD.code,1,3) %in% c("R95") , "Other ill-defined and unspecified causes of mortality",
        substr(ICD.code,1,1) %in% c("S") & as.numeric(substr(ICD.code,2,2)) <=1, "Injuries to the head and the neck",
        substr(ICD.code,1,2) %in% c("S2"), "Injuries to the thorax",
        substr(ICD.code,1,2) %in% c("S3"), "Injuries to the abdomen, lower back, lumbar spine and pelvis",
        substr(ICD.code,1,3) %in% c("S72") , "Fracture of femur",
        substr(ICD.code,1,1) %in% c("T") & as.numeric(substr(ICD.code,2,3)) >=20 & as.numeric(substr(ICD.code,2,3)) <=32, "Burns and corrosions",
        substr(ICD.code,1,4) %in% c("T391") , "Poisoning by 4-Aminophenol derivatives",
        substr(ICD.code,1,3) %in% c("T40") , "Poisoning by narcotics and psychodysleptics [hallucinogens,",
        substr(ICD.code,1,3) %in% c("T42") , "Poisoning by antiepileptic, sedative-hypnotic and antiparkinsonism drugs",
        substr(ICD.code,1,3) %in% c("T43") , "Poisoning by psychotropic drugs, not elsewhere classified",
        substr(ICD.code,1,4) %in% c("T509") , "Poisoning by other and unspecified drugs, medicaments and biological substances",
        substr(ICD.code,1,1) %in% c("T") & as.numeric(substr(ICD.code,2,3)) >=51 & as.numeric(substr(ICD.code,2,3)) <=65, "Toxic effects of substances chiefly nonmedicinal as to source",
        substr(ICD.code,1,3) %in% c("T58") , "Toxic effect of carbon monoxide",
        substr(ICD.code,1,3) %in% c("T71") , "Asphyxiation",
        substr(ICD.code,1,4) %in% c("T751") , "Drowning and nonfatal submersion",
        substr(ICD.code,1,4) %in% c("U071","U072") , "COVID-19",
        substr(ICD.code,1,1) %in% c("V") & as.numeric(substr(ICD.code,2,3)) >=1 & as.numeric(substr(ICD.code,2,3)) <=89, "Land transport accidents",
        substr(ICD.code,1,1) %in% c("W") & as.numeric(substr(ICD.code,2,3)) >=0 & as.numeric(substr(ICD.code,2,3)) <=19, "Falls",
        substr(ICD.code,1,1) %in% c("W") & as.numeric(substr(ICD.code,2,3)) >=65 & as.numeric(substr(ICD.code,2,3)) <=74, "Accidental drowning and submersion",
        substr(ICD.code,1,1) %in% c("X") & as.numeric(substr(ICD.code,2,3)) >=0 & as.numeric(substr(ICD.code,2,3)) <=9, "	Exposure to smoke, fire and flames",
        substr(ICD.code,1,3) %in% c("X41") , "Accidental poisoning by and exposure to antiepileptic, sedative-hypnotic, antiparkinsonism and psychotropic drugs, not elsewhere classified",
        substr(ICD.code,1,3) %in% c("X42") , "Accidental poisoning by and exposure to narcotics and psychodysleptics [hallucinogens,, not elsewhere classified",
        substr(ICD.code,1,3) %in% c("X44") , "Accidental poisoning by and exposure to other and unspecified drugs, medicaments and biological substances",
        substr(ICD.code,1,3) %in% c("X59") , "Accidental exposure to unspecified factor",
        substr(ICD.code,1,1) %in% c("X") & as.numeric(substr(ICD.code,2,3)) >=60 & as.numeric(substr(ICD.code,2,3)) <=84, "Intentional self-harm",
        substr(ICD.code,1,1) %in% c("X") & as.numeric(substr(ICD.code,2,3)) >=85, "Assault",
        substr(ICD.code,1,1) %in% c("Y") & as.numeric(substr(ICD.code,2,3)) <=9, "Assault",
        substr(ICD.code,1,1) %in% c("Y") & as.numeric(substr(ICD.code,2,3)) >=10 & as.numeric(substr(ICD.code,2,3)) <=34 , "Event of undetermined intent",
        !is.na(ICD.code), "Other")
}  

