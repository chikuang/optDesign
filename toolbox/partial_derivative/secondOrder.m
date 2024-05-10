function f = secondOrder(xi)
  f = [1, xi(1), xi(2), xi(1)^2, xi(2)^2, xi(1)*xi(2)]';
end
