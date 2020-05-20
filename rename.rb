puts "hi!"
puts "what file do you want to transform from edge to adjlist style?"
x = gets.chomp
lines = IO.readlines(x)

adj_list = {}
lines.each do |line|
    ary = line.split(' ').map(&:to_i)
    adj_list[ary[0]] ||= []
    adj_list[ary[1]] ||= []
    adj_list[ary[0]] += [ary[1]] 
end

ary = []
new_lines = adj_list.sort_by {|k,v| k}.map do |k, v|
    # k v0 v 1 v2
    t = [k, v.sort].flatten
end
puts 'life is good in ruby... '
fout = gets.chomp
f = File.new(fout, 'w')
IO.write(f, new_lines.map {|i| i.join(' ')}.join("\n"))


