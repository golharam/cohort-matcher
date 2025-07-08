cwlVersion: v1.2
class: ExpressionTool

inputs:
  bamsheet:
    type: File

outputs:
  sample_list:
    type:
      type: array
      items:
        type: record
        fields:
          sample_id: string
          bamfile: File

expression: >
  ${
    const fs = require("fs");
    const lines = fs.readFileSync(inputs.bamsheet.path, 'utf8').trim().split('\n').slice(1);
    return {
      sample_list: lines.map(line => {
        const [sample_id, bamfile_path] = line.split('\t');
        return {
          sample_id,
          bamfile: {"class": "File", "path": bamfile_path}
        };
      })
    };
  }
