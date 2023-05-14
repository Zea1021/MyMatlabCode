function index = RouletteWheelSelection(P)

    r = rand;
    
    index = find(r <= P, 1);

end