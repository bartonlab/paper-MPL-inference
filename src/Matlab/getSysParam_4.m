% this file contains the system parameters for various simulated datasets
if(strcmp(thisSet, 'medium_simple'))
	Lin = 50;
	DAll = 10;
    DeleAll = 10;
    selVal = 50;
    delSelVal = -50;
    Nin = 1000;
	T = 1000;
	Tstart = 0;
	Tused = 1000;
    dT = 10;
    ng = 100;
	muVal = 0.0001;
	perSiteSelction = [selVal/2/Nin*ones(1, DAll) zeros(1, Lin - (DAll+DeleAll)) delSelVal/2/Nin*ones(1, DeleAll)];
	classesOfSites = 3;
elseif(strcmp(thisSet, 'medium_complex'))
	Lin = 50;
	DAll = 10;
    DeleAll = 10;
    selVal = 200;
    delSelVal = -200;
    Nin = 1000;
	T = 1000;
	Tstart = 10;
	Tused = 300;
    dT = 10;
    ng = 100;
	muVal = 0.0001;
	perSiteSelction = [selVal/2/Nin*ones(1, DAll) zeros(1, Lin - (DAll+DeleAll)) delSelVal/2/Nin*ones(1, DeleAll)];
	classesOfSites = 3;
end
