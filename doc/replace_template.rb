# Replace the template with markdown contents
# @Author Xin Feng
# @Date 11/28/2014


if ARGV.length < 3
  $stderr.puts "Usage: ruby #{$0} template.html content.md target.html"
  exit
end
template = File.open(ARGV[0], "rb")
content = File.open(ARGV[1], "rb")
content = content.read
template = template.read
template.gsub!("GGPLAY",content)
File.open(ARGV[2], 'w') {|f| f.write(template) }
