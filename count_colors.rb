puts "hi"

lines = IO.readlines('colorings.txt').map(&:to_i)

keys = lines.uniq
values = Array.new(keys.size, 0)

h = Hash[keys.zip values]

lines.each do |x|
    h[x] +=1
end

vals = h.values.sort.reverse
puts vals

puts 'life is good in ruby... '
fout = gets.chomp
f = File.new(fout, 'w')
IO.write(f, vals[0..4].map(&:to_s).join(","))