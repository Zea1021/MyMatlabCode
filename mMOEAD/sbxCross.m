function [indiv1sub,indiv2sub] = sbxCross(indiv1, indiv2)
L = length(indiv1);
indiv1sub = zeros(1, L);
indiv2sub = zeros(1, L);
Rand = rand(1, L);
indiv1sub(Rand > 0.5)   = indiv1(Rand > 0.5);
indiv1sub(Rand <= 0.5)  = indiv2(Rand <= 0.5);
indiv2sub(Rand >= 0.5)  = indiv2(Rand >= 0.5);
indiv2sub(Rand < 0.5)   = indiv1(Rand < 0.5);

