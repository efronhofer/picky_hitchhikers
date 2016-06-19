{
									
	Simulations of Insects in Structured Populations			
			-----------------------------			
		 Field Station Fabrikschleichach		
			Evolutionary Ecology Group			
			  University of Wuerzburg			
									
	    	 

    This program is free software; you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation; either version 3 of the License, or
    (at your option) any later version.
        
    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.
        
    You should have received a copy of the GNU General Public License
    along with this program; if not, write to the Free Software
    Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston,
    MA 02110-1301, USA.
}
program SISP_spadiseius;

uses Math;

//------------------------------------------------------------------------------
//constants --------------------------------------------------------------------
const
	XMAX = 32;      													// no. of patches in the world (in x and y dim.)
	YMAX = 32;
	PH = 0.5;															// suitable haitat
	ITMAX = 150;														// no of iterations per patch for semicontinuous model
	KMAX = 1000;    													// maximal habitat capacity
	VMAX = 100;															// maximal vector capacity
	RS = 2;           													// random seed
	VARI = 0.2;       													// variance used for mutation
	PI = 3.141593;														// pi
	MAXRUNS = 10;														// no of repeats

//------------------------------------------------------------------------------
//types ------------------------------------------------------------------------
type
// One individual --------------------------------------------------------------
TFemale = object
	locus_pref : array[1..4] of single;                      			// locus for vector preference: Artibeus, Chasmodia, Cholus, Trigona
  end;

// One Patch -------------------------------------------------------------------
TPatch = object
	females  : array[1..KMAX] of TFemale;                               // all females
	no_females : word;                                               	// counts all females
	newFemales : array[1..KMAX] of TFemale;    		                    // all new females (newly dispersed or born)
	no_newFemales : word;                                            	// counts all new females
	suitable : byte;													// patch availabe or not
	age : byte;															// records patch age; info for second dispersal step
  end;

// World -----------------------------------------------------------------------
TWorld = array[1..XMAX,1..YMAX] of TPatch;                              //whole metapopulation

//------------------------------------------------------------------------------
// variables -------------------------------------------------------------------
var
	metapop : TWorld;             										// the whole metapopulation
	capacity : word;       	    										// habitat capacity
	lambda_null : single;           									// basic fertility i.e. mean no. of offspring
	sigma : single;                 									// environmental fluctuations
	maxPage : byte;          											// max. patch age
	mutationRate : single;          									// mutation rate
	generations : word;        											// no. of generations
	time : word;               											// counter for generations
	input, map, data_I, data_SSP : text; 								// data output and input files  
	run_no : word;														// counter for repeats
	metapopsize : longint;												// saves metapop size
	outputfilename_SSP : string;										// outfile
	dispCapacity : array[1..4] of word;									// capacity of vectors
	dispMeanDist : array[1..4] of word;									// mean dispersal distance of vectors
	dispMort : array[1..4] of single;									// mean dispersal mortality of vectors
	dispRate : array[1..4] of single;									// vector specific emigration rate
	dispPref : array[1..4] of string;									// vector specific preference for young patches
	help_probs : array[1..4] of real;									// probability of event, i.e. visits
	real_emi : longint;													// counters for emigrants

//------------------------------------------------------------------------------
// functions -------------------------------------------------------------------

// poisson distribution -----------------------------------------------------
function Poisson(lamb:double):integer;
// this function creates Poisson distributed random numbers
// with mean lamb. the source code can be found in:
// Press WH, Flannery BP, Teukolsky SA, Vetterling WT: 
// Numerical Recipes in Pascal. The Art of Scientific Computing. 
// Cambridge, New York, Melbourne: Cambridge University Press; 1989. 
end;

// log Gauss distribution ---------------------------------------------------
function loggauss(Fertility: real): real;
// calculates the mean random fertility
// of a single female individual from a
// Gaussian probability distribution,
// using the variable "Fertility" as
// mean value and "sigma" (environmental
// fluctuations) as standard deviation. the source code can be found in:
// Press WH, Flannery BP, Teukolsky SA, Vetterling WT: 
// Numerical Recipes in Pascal. The Art of Scientific Computing. 
// Cambridge, New York, Melbourne: Cambridge University Press; 1989. 
end;

// dispersal kernel ---------------------------------------------------------
function skewKernel(meanNE:integer):real;
//skewed dispersal kernel following 
// Gros A, Poethke HJ & Hovestadt T (2006) 
// Evolution of local adaptation in dispersal strategies. 
// Oikos 114: 544-552.
end;

// mutation -----------------------------------------------------------------
function mutate(const allele : real) : real;
	// introduces a mutation uniformely distributet around the mean "allele" with
	// variance VARI (see constants)!
	// from A. Poethke
	// called by procedure "Reproduction"
	var
		m1 : real;      
	begin
	m1 := allele + VARI * (2*random - 1);
	// limit mutation between 0 and 1
	if m1 < 0 then m1 := abs(m1);
	if m1 > 1 then m1 := 2 - m1;
	mutate := m1;
	end;

// offsping -------------------------------------------------------------------
function NoOffspring(const Nt: longint; lambdaLoc : real): integer;
	// calculates the actual no of offspring taking into account
	// logistic growth
	// called by procedure "reproduction"
	// Nt : actual pop size
	var
		survival,a : real;
	begin
	a:=(lambda_null-1)/capacity;											//calculation of crowding-parameter
	if a<0 then a:=0;	
	survival := 1/(1+a*Nt);
	NoOffspring:=poisson(lambdaLoc*survival);
	end;

// -----------------------------------------------------------------------------
// biol. procedures ------------------------------------------------------------

// informed patch choice ---------------------------------------------------------
procedure PatchChoice (const x,y : integer; const info : string; var xn, yn : integer);
	// procedure to find nearest neighbours
	procedure findNNX (const x,y,nnd_no: integer; var xn, yn : integer);
		begin
			case nnd_no of   //3
			1 : begin
				xn := x - 1;
				yn := y + 1;
				end;
			2:  begin
				xn := x;
				yn := y +1;
				end;
			3:  begin
				xn := x +1;
				yn := y +1;
				end;
			4:  begin
				xn := x -1;
				yn := y;
				end;
			5:  begin
				xn := x +1;
				yn := y;
				end;
			6:  begin
				xn := x -1;
				yn := y-1;
				end;
			7:  begin
				xn := x;
				yn := y-1;
				end;
			8:  begin
				xn := x+1;
				yn := y-1;
				end;
			end;   //3

		if xn < 1 then xn := XMAX;
		if xn > XMAX then xn := 1;
		if yn < 1 then yn := YMAX;
		if yn > YMAX then yn := 1;
		end;
	
var
	nnd_no,n : integer;
	age_nnb4, min_nn : integer; 
	posx, posy : integer;
	noInfoPos : array[1..8] of integer;
	length_noInfoPos : integer;
   
begin   //1
	if info = 'true' then
	begin//2
		for n:= 1 to 8 do
		begin//3
			// determine possible new patch
			findNNX(x,y,n,posx,posy);
			if n = 1 then 
			begin //4
				age_nnb4 := metapop[posx,posy].age;
				min_nn := 1;
			end //4
			else
			begin//4
				if metapop[posx,posy].age < age_nnb4 then
				begin //5
					age_nnb4 := metapop[posx,posy].age;
					min_nn := n
				end; //5
			end;//4
			
		end;//3
	// go to youngest patch
	findNNX(x,y,min_nn,xn,yn);
	end//2
	else
	// new patch address if no information is used
	begin //2
		length_noInfoPos := 0;
		// for comparability the vector still prefers suitable over non-suitable
		for n:= 1 to 8 do
		begin//3
			// determine possible new patch
			findNNX(x,y,n,posx,posy);
			// check for suitability
			if (metapop[posx,posy].suitable=1) then
			begin //4
				inc(length_noInfoPos);
				noInfoPos[length_noInfoPos] := n;
			end; //4
		end; //3	
	
		if (length_noInfoPos > 0) then nnd_no := noInfoPos[random(length_noInfoPos) + 1]
		else nnd_no := random(8)+1;
		// determine new patch
		findNNX(x,y,nnd_no,xn,yn);
	end; //2
end;  //1

// dispersal ----------------------------------------------------------------
procedure dispersal;
	var
	patchNo_x, patchNo_y : integer;
	e_help : real;
	its,f, help, strat, d : integer; 									
	cap_cnt, focal : integer;
	dispDist, phi : real;
	xstart, ystart : real;
	xend, yend : integer;
	x2end, y2end : integer;
	thisgroup : array[1..VMAX] of TFemale;								
	help_fs : array[1..KMAX] of integer;							
	help_cnt : integer;
	help_focal,i : integer;


	begin    //1
	// reset realized emigrants counter
	real_emi:=0;
	
	for patchNo_x := 1 to XMAX do             //loops for every patch
	for patchNo_y := 1 to YMAX do
	begin //2
	
		// check whether inflorescence is suitable and inhabited
		if ((metapop[patchNo_x,patchNo_y].suitable=1) and (metapop[patchNo_x,patchNo_y].no_females > 0)) then
		begin //3
			
			// add loop for iterations
			for its:= 1 to ITMAX do
			begin // 4
				
				//check pop size
				if (metapop[patchNo_x,patchNo_y].no_females < 1) then break;
				
				// here: randomly choose an event (strat)
				// possibilities are: (1) Artibeus, (2) Chasmodia, (3) Cholus, (4) Trigona, or nothing happens (0)
				strat:=0;
				e_help := random;
				
				// check prefernces
				if (e_help < help_probs[1]) then 																				// Artibeus
				begin // 5
					strat:=1;
				end; // 5																																				
				if ((e_help >= help_probs[1]) and (e_help < help_probs[2])) then		// Chasmodia
				begin // 5
					strat:=2;
				end; // 5
				if ((e_help >= help_probs[2]) and (e_help < help_probs[3])) then 		// Cholus
				begin // 5
					strat:=3;
				end; // 5
				if ((e_help >= help_probs[3]) and (e_help < help_probs[4])) then 		// Trigona
				begin // 5
					strat:=4;
				end; // 5
			
				// now that the event is chosen:
				if not(strat = 0) then
				begin //5
				
				// help vector for animal adresses;
				for i:=1 to metapop[patchNo_x,patchNo_y].no_females do help_fs[i] := i;
				help_cnt := metapop[patchNo_x,patchNo_y].no_females;
				
					// reset vector capacity counter
					cap_cnt := 0;
					// here comes a loop over all x indivs that fit on the vector "strat" or until the patch is empty
					repeat //6
						if help_cnt < 1 then break;
						// randomly choose an indiv from patch
						help_focal := random(help_cnt)+1;
						focal := help_fs[help_focal];
						
						help_fs[help_focal] := help_fs[help_cnt];
						dec(help_cnt);

							if (random < metapop[patchNo_x,patchNo_y].females[focal].locus_pref[strat]) then
							begin // 8
								// inc. counter for vect capacity
								inc(cap_cnt);
								// add to actual group
								thisgroup[cap_cnt] := metapop[patchNo_x,patchNo_y].females[focal];
								// delete indiv. from patch
								metapop[patchNo_x,patchNo_y].females[focal] := metapop[patchNo_x,patchNo_y].females[metapop[patchNo_x,patchNo_y].no_females];
								// dec counter for remaning indivs in patch
								dec(metapop[patchNo_x,patchNo_y].no_females);
								// counter for emigraion rate
								inc(real_emi);
							end; // does the animal take the focal vector
			
					until ((cap_cnt = dispCapacity[strat])or(help_cnt < 1));		// until group is full or patch is empty
					
					// first calculate distance
					dispDist := skewKernel(dispMeanDist[strat]);

					//random direction
					phi := 2 * PI * random;    

					//area-to-area dispersal
					xstart := patchNo_x + random - 0.5;    
					ystart := patchNo_y + random - 0.5;

					// target coordinates
					xend := round( xstart + dispDist * cos(phi) );    
					yend := round( ystart + dispDist * sin(phi) );
									
					// make torus
					if xend > XMAX then
					begin//6
						xend := xend mod XMAX;
						if xend = 0 then xend := XMAX;
					end;//6

					if xend < 1 then
					begin//6
						xend := XMAX - (abs(xend) mod XMAX);
					end;//6

					if yend > YMAX then
					begin//6
						yend := yend mod YMAX;
						if yend = 0 then yend := YMAX;
					end;//6

					if yend < 1 then
					begin//6
						yend := YMAX - (abs(yend) mod YMAX);
					end;//6
							
					// here include a second dispersal step: NN8 and use of information on patch age
					PatchChoice(xend, yend, dispPref[strat], x2end, y2end);
						
					// now that target is clear do dispersal and specific mortality
					// dispersal if mortality does not apply and patch is suitable; copy indivs. into new indivs vector of new patch!							
					if ((random > dispMort[strat]) and (metapop[x2end,y2end].suitable=1)) then
					begin //7
						for d:=1 to cap_cnt do
						begin //6					
							// put dispersers into new patch
							inc(metapop[x2end,y2end].no_newFemales);	
							metapop[x2end,y2end].newFemales[metapop[x2end,y2end].no_newFemales] := thisgroup[d];
						end; //7
					end; //6
				end; //5 if strat is not zero
			end; //4 end of loop for iterations!
		end; //3
	end; //2 patch loop
	
	//now that dispersal is over for all patches, merge dispersers (newFemales) and residents (females)
	for patchNo_x := 1 to XMAX do           //loops for every patch
	for patchNo_y := 1 to YMAX do
	begin   //2
		// if there are females and patch is suitable...
		if ((metapop[patchNo_x,patchNo_y].no_newFemales > 0) and (metapop[patchNo_x,patchNo_y].suitable=1)) then
		begin //3
			help := metapop[patchNo_x,patchNo_y].no_females; // counter for array newFemales
			for f:= 1 to metapop[patchNo_x,patchNo_y].no_newFemales do      // loop for every female
			begin //4
				inc(help);
				// insert new  females in resident female array
				metapop[patchNo_x,patchNo_y].females[help] := metapop[patchNo_x,patchNo_y].newFemales[f];
			end; //4
			// increase counter for no of females
			metapop[patchNo_x,patchNo_y].no_females := help;
			//reset counter for new indiv
			metapop[patchNo_x,patchNo_y].no_newFemales := 0;
		end; //3
	end; //2
end; //1

// reproduction -------------------------------------------------------------
procedure reproduction;
	var
		patchNo_x, patchNo_y : integer;        // coordinates of actual patch
		f,o,p : integer;                       // counter for females, males and offspring
		no_offspring : integer;                // no. of offspring
		pop_size : longint;                    // actual pop. size
		lambdaPatch : real;

	begin //1
		for patchNo_x := 1 to XMAX do             //loops for every patch
		for patchNo_y := 1 to YMAX do
		begin   //2
			//reset counter for new indiv
			metapop[patchNo_x,patchNo_y].no_newFemales := 0;
			//check whether there are females
			if ((metapop[patchNo_x,patchNo_y].suitable=1) and (metapop[patchNo_x,patchNo_y].no_females > 0)) then
			begin //3
				// calculate actual pop. size
				pop_size := metapop[patchNo_x,patchNo_y].no_females;
				// calculate lambda after incorporation of var sigma (environmental fluctuation)
				if lambda_null > 0 then lambdaPatch:= loggauss(lambda_null)  
				else lambdaPatch := 0;
				//loop for every female
				for f:=1 to metapop[patchNo_x,patchNo_y].no_females do
				begin //4
					// calculate no. of offspring
					// this function incorporates logistic growth
					// then the no. of offspring is calculated from a poisson distribution
					no_offspring := NoOffspring(pop_size,lambdaPatch);
					if no_offspring > 0 then
					begin //5
						for o := 1 to no_offspring do
						begin //6
								inc(metapop[patchNo_x,patchNo_y].no_newFemales);
								// determine genetics of offspring
								// mutation rate is introduced
								// vector preference locus
								for p:= 1 to 4 do
								begin //7
								metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_pref[p] := metapop[patchNo_x,patchNo_y].females[f].locus_pref[p];
								if random < mutationRate then metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_pref[p] := mutate(metapop[patchNo_x,patchNo_y].newFemales[metapop[patchNo_x,patchNo_y].no_newFemales].locus_pref[p]);
								end; //7
						end; //6
					end;  //5
				end; //4
			end;  //3
		end; //2
	end; //1

// death of organisms -------------------------------------------------
procedure death;
var
	patchNo_x, patchNo_y : integer;     // coordinates of actual patch
	f 					: integer; 		// females counter
begin // 1
	// reset metapopsize
	metapopsize := 0;
	for patchNo_x := 1 to XMAX do   // loop in x dir.
    for patchNo_y := 1 to YMAX do  // loop in y dir.
    begin //2
		if (metapop[patchNo_x,patchNo_y].suitable = 1) then
		begin //3
		// replace adult female vector with offspring vector
				for f:= 1 to metapop[patchNo_x,patchNo_y].no_newFemales do
				begin //4
					metapop[patchNo_x,patchNo_y].females[f] := metapop[patchNo_x,patchNo_y].newFemales[f];
				end;  //4
				//reset female counters
				metapop[patchNo_x,patchNo_y].no_females := metapop[patchNo_x,patchNo_y].no_newFemales;
				metapop[patchNo_x,patchNo_y].no_newFemales := 0;
				
				// age patches
				inc(metapop[patchNo_x,patchNo_y].age);
				// deterministic patch extinction
				if metapop[patchNo_x,patchNo_y].age >= maxPage then
				begin //4
					metapop[patchNo_x,patchNo_y].no_females := 0;
					metapop[patchNo_x,patchNo_y].age := 0;
				end;  //4
				metapopsize := metapopsize + metapop[patchNo_x,patchNo_y].no_females;
		end; //3
	end; // 2
	
	// some console output
	if (time mod 10) = 0 then 
	begin
		write(run_no:10); write(time:10); write(metapopsize:10); writeln(real_emi:10);
	end;

end; // 1

// results output ---------------------------------------------------------------
procedure resultsOutput;
var
	patchNo_x, patchNo_y, f,p, occp, suitp : integer;     // coordinates of actual patch
	outputfilename_I,  runstring : string;

begin //1

	//individual output
    str(run_no, runstring);
    outputfilename_I := './results/output_Individuals_repeat' + runstring + '.out';
    assign(data_I, outputfilename_I);
    rewrite(data_I);

	writeln(data_I,'x      y       suitable	       Artibeus        Chasmodia       Cholus          Trigona');

	// reset occupancy
	occp := 0;
	suitp := 0;
	
	for patchNo_x := 1 to XMAX do   // loop in x dir.
    for patchNo_y := 1 to YMAX do  // loop in y dir.
    begin //2
    
		if metapop[patchNo_x,patchNo_y].no_females > 0 then inc(occp);
		if metapop[patchNo_x,patchNo_y].suitable = 1 then inc(suitp);    
		for f:= 1 to metapop[patchNo_x,patchNo_y].no_females do
		begin //3
			
			write(data_I,patchNo_x:10, '    ',patchNo_y:10, '    ',metapop[patchNo_x,patchNo_y].suitable:10, '       ');
			for p:= 1 to 4 do
			begin //4
				if p < 4 then write(data_I, metapop[patchNo_x,patchNo_y].females[f].locus_pref[p]:5:10, '    ')
				else writeln(data_I, metapop[patchNo_x,patchNo_y].females[f].locus_pref[p]:5:10);
			end; //4
		
		end; //3
    
    end; //2
	close(data_I);
    
    // SSP level output
    if run_no = 1 then writeln(data_SSP, 'run_no		time		sspsize			suitable		occupancy		emino');
    
    write(data_SSP, run_no:6, time:10, metapopsize:10, suitp:10, occp:10);
	writeln(data_SSP,'       ',real_emi:10);

end; //1

//--------------------------------------------------------------------------------
// initialize --------------------------------------------------------------------
procedure Initialize;
	var
		patchNo_x, patchNo_y : integer;
		f,p: integer;
		mapname, sizestring, phstring : string;
		i,j : integer;

	begin //1      
	
		// read parameters
		assign(input, 'input_SISP_V0.in');
		reset(input);
        
        readln(input); readln(input); readln(input); readln(input);    
		readln(input, capacity); readln(input);
		readln(input, lambda_null); readln(input);
		readln(input, sigma); readln(input);
		readln(input, maxPage); readln(input);	
		readln(input, mutationRate); readln(input);
		readln(input, generations); readln(input); readln(input);
		readln(input, dispCapacity[1]);readln(input);
		readln(input, dispCapacity[2]);readln(input);
		readln(input, dispCapacity[3]);readln(input);
		readln(input, dispCapacity[4]);readln(input);
		readln(input, dispMeanDist[1]);readln(input);
		readln(input, dispMeanDist[2]);readln(input);
		readln(input, dispMeanDist[3]);readln(input);
		readln(input, dispMeanDist[4]);readln(input);
		readln(input, dispMort[1]);readln(input);
		readln(input, dispMort[2]);readln(input);
		readln(input, dispMort[3]);readln(input);
		readln(input, dispMort[4]);readln(input);
		readln(input, dispRate[1]);readln(input);
		readln(input, dispRate[2]);readln(input);
		readln(input, dispRate[3]);readln(input);
		readln(input, dispRate[4]);readln(input);
		readln(input, dispPref[1]);readln(input);
		readln(input, dispPref[2]);readln(input);
		readln(input, dispPref[3]);readln(input);
		readln(input, dispPref[4]);
		
		close(input);
		
		// read landscape
		str(XMAX, sizestring);
		str(PH:2:1, phstring);
		mapname := './map/map' + sizestring + '_' + phstring + '.in';
		assign(map, mapname);
		reset(map);
		
		// calc probabilities for event occurrance
		help_probs[1] :=  dispRate[1];
		help_probs[2] :=  dispRate[1] + dispRate[2];
		help_probs[3] :=  dispRate[1] + dispRate[2] + dispRate[3];
		help_probs[4] :=  dispRate[1] + dispRate[2] + dispRate[3] + dispRate[4];
		
		//initialize patches in world
        for patchNo_x := 1 to XMAX do   // loop in x dir.
        for patchNo_y := 1 to YMAX do  // loop in y dir.
        begin   //2
			// initialize all patches empty
			metapop[patchNo_x,patchNo_y].no_females := 0;
			metapop[patchNo_x,patchNo_y].no_newFemales := 0;
			
			readln(map,metapop[patchNo_x,patchNo_y].suitable);
			
			if metapop[patchNo_x,patchNo_y].suitable = 1 then metapop[patchNo_x,patchNo_y].age := random(maxPage)
			else metapop[patchNo_x,patchNo_y].age := 99;

			//initialize females
			if metapop[patchNo_x,patchNo_y].suitable = 1 then metapop[patchNo_x,patchNo_y].no_females := capacity;

			for f := 1 to metapop[patchNo_x,patchNo_y].no_females do
			begin //3
				// random initialization of dispersal prefernce locus
				for p:=1 to 4 do
				begin
					metapop[patchNo_x,patchNo_y].females[f].locus_pref[p] := random;
				end;
			end; //3
		end;  //2
		
		close(map);
		
	end;  //1

// simulation ---------------------------------------------------------------------
// main ---------------------------------------------------------------------------
begin //1
	// intialize random number generator
	randseed := RS;  
	// output SSP level
    outputfilename_SSP := './results/output_SSP.out';
    assign(data_SSP, outputfilename_SSP);
    rewrite(data_SSP);

	for run_no:=1 to MAXRUNS do
	begin //2    
		Initialize;														// initializes world and reads form input file
		
		for time:= 1 to generations do             						// time loop
		begin //3
			dispersal;                          						// juvenile dispersal
			reproduction;                       						// reproduction in new patch
			death;                              						// death of adults and reset all dispersal makers for survivors
			if time = generations then resultsOutput;
		end; //3
	
	end; //2
	close(data_SSP);
end. //1
