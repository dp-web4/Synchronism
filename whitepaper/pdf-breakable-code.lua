-- Long inline code in the PDF gets break points.
--
-- pandoc sets inline code as \texttt{...}, which TeX never breaks: the mono
-- font has no hyphen character and a path has no spaces. The provenance this
-- paper cites inline (Research/proposals/<long_snake_name>.md,
-- synchronism-site/explorer/findings/<long-kebab-name>, scripts) therefore ran
-- straight off the page. Measured 2026-10-31: 98 words on 55 of 293 pages
-- ended past the 1in margin, most of them clipped at the physical page edge,
-- so the reader could not see the end of the path.
--
-- For inline code of MIN_LEN characters or more, this emits the \texttt
-- itself and puts \allowbreak after each / _ - . so a line can end there. Nothing
-- is inserted into the text, so the text layer is unchanged. Short code
-- (1/Ω_m, B_max) is left to pandoc, since a break there would split a symbol.
-- PDF only: make-pdf.sh passes this filter, and the HTML builds do not.

local MIN_LEN = 20

local escapes = {
  ['\\'] = '\\textbackslash{}',
  ['{'] = '\\{', ['}'] = '\\}',
  ['#'] = '\\#', ['$'] = '\\$', ['%'] = '\\%', ['&'] = '\\&', ['_'] = '\\_',
  ['~'] = '\\textasciitilde{}', ['^'] = '\\^{}',
  ['`'] = '\\textasciigrave{}', ["'"] = '\\textquotesingle{}',
  [' '] = '\\ ',
}
local breakafter = { ['/'] = true, ['_'] = true, ['-'] = true, ['.'] = true }

function Code(el)
  if not FORMAT:match('latex') then return nil end
  if utf8.len(el.text) == nil or utf8.len(el.text) < MIN_LEN then return nil end
  local out = {}
  for _, cp in utf8.codes(el.text) do
    local ch = utf8.char(cp)
    out[#out + 1] = escapes[ch] or ch
    -- '-' gets '{}' first so a following '-' cannot form a ligature
    if ch == '-' then out[#out + 1] = '{}' end
    if breakafter[ch] then out[#out + 1] = '\\allowbreak{}' end
  end
  return pandoc.RawInline('latex', '\\texttt{' .. table.concat(out) .. '}')
end
