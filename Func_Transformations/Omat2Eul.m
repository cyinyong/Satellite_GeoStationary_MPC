function eul = Omat2Eul(Omat)

eul(2) = asin(-Omat(1,3));
eul(1) = asin(Omat(2,3)/cos(eul(2)));
eul(3) = asin(Omat(1,2)/cos(eul(2)));

end